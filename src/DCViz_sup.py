# -*- coding: utf-8 -*-

import re, numpy, time, sys, signal, os, struct
import matplotlib.pylab as plab
from os.path import join as pjoin


class DCVizPlotter:
    
    figMap = {}

    armaBin = False
    
    isFamilyMember = False
    familyName = "unnamed"
    familyFileNames = []

    skippedRows = []
    skipRows = None

    delay = 3

    parent = None
    
    
    canStart = False
    
    def __init__(self, filepath=None, dynamic=False, useGUI=False):
        self.dynamic = dynamic
        self.useGUI = useGUI
        
        self.plotted = False
        self.stopped = False
        self.SIGINT_CAPTURED = False
            
        self.figures = []
    
        self.filepath = filepath
        self.file = None
        
        signal.signal(signal.SIGINT, self.signal_handler)
    
    
    def signal_handler(self, signal, frame):
        print "Ending session..."
        self.SIGINT_CAPTURED = True
    
    def __str__(self):
        
        if self.isFamilyMember:
            return self.familyName + " (family)"
        return ".".join(os.path.split(self.filepath)[-1].split(".")[0:-1])

        
    def get_data(self, setUpFamily):
    
        if setUpFamily:
            
            familyHome, myName = os.path.split(self.filepath)
            familyNames = [name for name in os.listdir(familyHome)\
                            if name != myName and re.findall(self.nametag, name)]
            
            familyMembers = [os.path.join(familyHome, name) for name in familyNames]
            
            data = [0]*(len(familyMembers) + 1)
            self.familyFileNames = [0]*len(data)     
            
            self.familyFileNames[0] = myName
            data[0] = self.get_data(setUpFamily=False)
            
            for i in range(len(familyMembers)):
             
                self.file = None
                self.filepath = familyMembers[i]
                
                self.familyFileNames[i+1] = os.path.split(familyMembers[i])[1]
                data[i+1] = self.get_data(setUpFamily=False)
                
            self.file.close()
            return data
            
        self.reload()
        
        if self.armaBin:
            data = self.unpackArmaMatBin(self.file)            
        else:
            data = numpy.array(self.rx.findall(self.file.read()), numpy.float)

        self.file.close()
        
        output = [0]*data.shape[1]
        for i in range(data.shape[1]):
            output[i] = data[:,i]
        
        return tuple(output)
    
    def unpackArmaMatBin(self, armaFile):   

        unpacker = struct.Struct("d")
    
        armaFormat = armaFile.readline().strip()
        size = int(armaFormat[-1])
        
        n, m = armaFile.readline().strip().split()
        n = int(n)       
        m = int(m)
        
        data = numpy.zeros(shape=(n, m))
        
        for i in range(n):
            for j in range(m):                
                data[i][j] = unpacker.unpack(armaFile.read(size))[0]
            
        return data
    
    def set_figures(self):
        
        s = ""
        i = 0
        self.figures = []
        for fig in self.figMap.keys():
            s += "self.%s = plab.figure(%d); " % (fig, i)
            s += "self.i%s = self.add_figure(self.%s); " % (fig, fig)
        
            subFigs = self.figMap[fig]
            nFigs = len(subFigs)
            
            for j in range(nFigs):
                s += "self.%s = self.%s.add_subplot(%d, 1, %d); " % (subFigs[j], fig, nFigs, j+1)
                s += "self.add_subfigure(self.%s, self.i%s); " % (subFigs[j], fig)
            
            i += 1
        exec(s)
      
      
    def manageFigures(self):
        if not self.useGUI:
            self.set_figures()
        else:
            
            self.parent.comm.setFigureSignal.emit()
            
            i = 0
            while not self.canStart:
                time.sleep(0.01)
                i+=1
                if i > 500:
                    self.Error("TIMEOUT: Figures wasn't set...")
                    return          


    def Error(self, s):
        if self.useGUI:
            self.parent.parent.terminalTracker("DCViz", s)
        else:
            print s


    def waitForGreenLight(self):

        if self.dynamic:
            light = 'red'
            while light == 'red':
                
                if self.stopped or self.SIGINT_CAPTURED:
                    return True
                    
                light = self.load_sample()
                time.sleep(1)
        else:
            light = self.load_sample()
            if light == 'red':
                return True
        
        return False
                
    def shouldReplot(self):
        return not self.stopped and not self.SIGINT_CAPTURED
    def shouldBreak(self):
        return self.stopped or self.SIGINT_CAPTURED

    def showFigures(self):
        
        if not self.useGUI:       
            self.show()
        else:
            if self.plotted and self.dynamic:
                try:
                    self.show(drawOnly=True)
                except:
                    #Exception needed in order for 
                    #the thread to survive being closed by GUI
                    pass
            if not self.plotted:
                self.parent.comm.plotSignal.emit()

    def sleep(self):
        for i in range(int(self.delay)):
            time.sleep(1)
            if self.shouldBreak():
                break

        time.sleep(self.delay - int(self.delay))

    def mainloop(self):
        
        if not self.armaBin:
            breakMe = self.waitForGreenLight()
            
            if breakMe:
                return

        self.manageFigures()

        while (self.shouldReplot()):
       
            self.clear()
           
            data = self.get_data(setUpFamily = self.isFamilyMember)
          
            self.plot(data)  
            self.showFigures()
            self.plotted = True

                
            if self.dynamic:
                self.sleep()     
            else:
                if not self.useGUI:
                    raw_input("Press any key to exit")
                break
                
        if not self.useGUI:
            self.close()
        
            
            
    def load_sample(self):
        
        self.reload()
  
        sample = self.file.read()
        
        if not sample:
            self.Error("No data in file...")
            return "red"
            
        skipRows, self.N = self.sniffer(sample)

        self.file.seek(0)
        
        self.skippedRows = []
        for i in range(skipRows):
            self.skippedRows.append(self.file.readline().strip())
        
        anyNumber = r'[\+\-]?\d+\.?\d*[eE]?[\+\-]?\d*'
        self.rx = re.compile((r'(%s)\s+' % anyNumber)*(self.N-1) + r'(%s)[\n$]' % anyNumber)    
        
        return "green"
        
    
    def sniffer(self, sample):
        
        sampleList = [row.split() for row in sample.split("\n")]

        #If user has specified the number of rows to skip
        if self.skipRows is not None:
            print "skipped rows is given."
            return self.skipRows, len(sampleList[self.skipRows])

        nLast = len(sampleList[-1])
        i = len(sampleList) - 1

        if nLast != len(sampleList[-2]):
            nLast = len(sampleList[-2])
            i -= 1
        
        while len(sampleList[i]) == nLast and i >= 0:
            i -= 1

        #(nRows to skip = i+1, nCols = nLast)
        return i+1, nLast    
    
    def reload(self):
    
        if self.file:
            if not self.file.closed: 
                self.file.close()
            
        self.file = open(self.filepath , "r")

        if not self.armaBin and self.skipRows is not None:        
            for i in range(self.skipRows):
                self.file.readline()
        
        
        
        
            
    def add_figure(self, fig):
        self.figures.append([fig])
        return len(self.figures) - 1
        
    def add_subfigure(self, subfig, i):
        self.figures[i].append(subfig)
        return len(self.figures[i]) - 1
        
    def show(self, drawOnly=False):
        for fig in self.figures:
            fig[0].canvas.draw()
            if not drawOnly:
                fig[0].show()
        
    def clear(self):
        for fig in self.figures:
            for subfig in fig[1:]:
                subfig.clear()
                subfig.axes.clear()
                subfig.axes.legend_ = None
        
    def close(self):
        self.clear()
        
        for fig in self.figures:
            plab.close(fig[0])
        
        if self.file is not None:
            if not self.file.closed:
                self.file.close()
        
    def plot(self, data):
        "I am virtual"
    
    
