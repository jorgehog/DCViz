# -*- coding: utf-8 -*-

import re, numpy, time, signal, os, sys, struct, itertools

import matplotlib.pylab as plab

from matplotlib import rcParams
from os.path import join as pjoin
from random import shuffle

class dataGenerator:
    def __init__(self, data):
        self.data = data
        
        if len(data.shape) == 2:
            self.n, self.m = data.shape
            self.getD = self.get2Ddata
        elif len(data.shape) == 3:
            self.n, self.m, self.l = data.shape
            self.getD = self.get3Ddata
        else:
            self.m, = data.shape
            self.n = 1
            self.getD = self.get1Ddata
        
        self.shape = data.shape
        self.fullshape = (self.n, self.m)
        self.size = data.size

    def get1Ddata(self, i):
        if i == slice(0, 9223372036854775807, None):
            return self.data
        elif i == 0:
            return self.data
        elif i == slice(0, 1, None):
            return self.data
        else:
            raise IndexError("Index out of bounds.")
    
    def get3Ddata(self, i):
        raise NotImplementedError()

    def get2Ddata(self, i):
        return self.data[:, i]

    def __iter__(self):

        for i in range(self.m):
            yield self.getD(i)

    def __len__(self):
        return self.n

    def __getitem__(self, i):
        return self.getD(i)

    def __str__(self):
        return str(self.data)

    def __add__(self, other):
            
        if isinstance(other, dataGenerator):
            return self.data + other.data
        else:        
            raise NotImplementedError("add/sub only works for two datasets. Use 'sum([sum(d) for d in data])'")

    def __sub__(self, other):    
            return self.data + (-other.data)


class DCVizPlotter:
    
    figMap = {}
    
    nametag = None

    armaBin = False
    numpyBin = False    
    fileBin = False
    
    path = None
    filename = None
    filepath = None
    
    binaryHeaderBitSizes = None
    binaryHeader = None
    nColsFromHeaderLoc = None
    binaryHeaderTypes = None
    defaultBinaryHeaderTypes = {4: 'i', 8: 'd'}
    
    Ncols = None
    transpose = False    
    
    isFamilyMember = False
    familyName = None
    familyFileNames = []
    familyHome = None
    loadLatest = False
    loadSequential = False
    getNumberForSort = lambda _s, x : int(re.findall(_s.nametag, x)[0])
    ziggyMagicNumber = 1
    smartIncrement = True
    originalFilename = None


    skippedRows = []
    skippedCols = []
    skipCols = 0
    skipRows = None
    
    nFig = 0

    delay = 3

    parent = None
    
    stack = "V"
    
    canStart = False
    
    gifLoopDelay = 0
    
    transparent = False
 
    hugifyFonts = False
    labelSize= 20
    fontSize = 20 #Only invoked by hugifyFonts = True
    tickSize = 2
    
    fig_size = None
    relative_fig_size = None
    
    markers = itertools.cycle(['*', 'o', '^', '+', 'x'])
    lines   = itertools.cycle(['-', '--'])
    colors  = itertools.cycle(['b', 'r', 'g', 'k', 'c'])

    anyNumber = r'[\+\-]?\d+\.?\d*[eE]?[\+\-]?\d*|[\+\-]?nan|[\+\-]?inf'    

    
    def __init__(self, filepath=None, dynamic=False, useGUI=False, toFile=False, threaded=False):

        if not self.familyName:
            self.familyName = self.__class__.__name__

        if not self.nametag:
            raise RuntimeError("An instance of DCViz must have a specified nametag.")
        
        if self.fig_size:
            plab.rcParams['figure.figsize'] = self.fig_size
        elif self.relative_fig_size:
            plab.rcParams['figure.figsize'] = [old*scale for old, scale in zip(plab.rcParams['figure.figsize'], self.relative_fig_size)]
        
        self.dynamic = dynamic
        self.useGUI = useGUI
        
        self.plotted = False
        self.stopped = False
        self.SIGINT_CAPTURED = False
        self.makeGif = False
        self.end_seq = False
        
        self.figures = []
        
        self.filepath = filepath
        self.file = None
        
        self.toFile = toFile
        
        if self.loadSequential:
            self.nextInLine = 0
        
        if not (threaded or useGUI) and dynamic:
            print "[  DCViz  ]", "Interrupt dynamic mode with CTRL+C"
            signal.signal(signal.SIGINT, self.signal_handler)
        
        if self.Ncols is None and self.fileBin and not (self.binaryHeaderBitSizes and self.nColsFromHeaderLoc):
            self.Error("You need to specify the number of cols 'Ncols' in order to read a binary file.")
            self.Exit()
        
    
    
    def signal_handler(self, signal, frame):
        print "[%s] Ending session..." % "DCViz".center(10)
        self.SIGINT_CAPTURED = True
    
    def __str__(self):
        
        if self.isFamilyMember:
            return self.familyName + " (family)"
        return ".".join(os.path.split(self.filepath)[-1].split(".")[0:-1])

    def getRandomStyles(self):
        styles =  [col + line + mark for col in self.colors for line in self.lines for mark in self.markers]
        shuffle(styles)
        return styles

    def get_family(self):
        
        familyHome, _ = os.path.split(self.filepath)
            
        if not familyHome:
            self.Error("Single file name given as path. DCViz needs the absolute file path to work with families.")
            self.Exit()
        
        familyNames = [name for name in os.listdir(familyHome)\
                        if re.findall(self.nametag, name) and os.path.exists(pjoin(familyHome, name)) and "tmp" not in name]
       
        familyMembers = sorted([pjoin(familyHome, name) for name in familyNames])     
        
        return familyMembers

    def dictify(self, data):

        _data = {}        
        for fdata, name in zip(data, self.familyFileNames):
            _data[name] = fdata
            
        data = _data

    def get_data(self, setUpFamily):

        if setUpFamily:

            familyMembers = self.get_family()
            
            N = len(familyMembers)            

            if self.loadLatest:
                
                _file = familyMembers[0]
                oldest = os.path.getctime(_file)
                
                for member in familyMembers:
                    
                    try:
                        age = os.path.getctime(member) 
                        if  age > oldest:
                            oldest = age
                            _file = member
                    except:
                        pass
                        
                        
                self.filepath = _file
                self.file = None
                
                return self.get_data(setUpFamily=False)                
                
                
            elif self.loadSequential:
                
                
                if self.getNumberForSort:
                    
                    familyMembers = sorted(familyMembers, key=self.getNumberForSort)
                    
                    
                self.filepath = familyMembers[self.nextInLine]
      
                prev = self.nextInLine
                
                
                self.nextInLine = (self.nextInLine + self.ziggyMagicNumber)%N
                
                    
                
                if (self.nextInLine < prev and self.makeGif) or self.end_seq:
                    
                    if self.end_seq:
                        self.dynamic = False
                                        
                    self.end_seq = True
                    
                    if self.toFile:
                        print "WARNING: ToFile and makeGif should not be used together."
                
                self.file = None
                
                return self.get_data(setUpFamily=False)                
                
                
            else:
                
                self.familySkippedRows = []

                if self.fileBin and self.binaryHeaderBitSizes:
                    self.familyHeaders = []

                data = [0]*N
                self.familyFileNames = [0]*len(data)     
                
                for i in range(N):
                 
                    self.file = None
                    self.filepath = familyMembers[i]
                
                    self.familyFileNames[i] = os.path.basename(familyMembers[i])
                    
                    data[i] = self.get_data(setUpFamily=False)
                    
                    if self.skipRows:
                        self.familySkippedRows.append(self.skippedRows)
                    
                    if self.fileBin and self.binaryHeaderBitSizes:
                        self.familyHeaders.append(self.binaryHeader)

    

                    
            self.file.close()
            return data
            
        data = []
        #attempt to reload file untill data is found.
        t0 = time.time()
        while len(data) == 0:

            self.reload()

            if self.armaBin:
                data = self.unpackArmaMatBin(self.file)
            elif self.fileBin:
                data = self.unpackBinFile(self.file)
            elif self.numpyBin:
                data = self.unpackNumpyBin(self.file)
            else:
                data = []

                for line in self.file:
                    split = line.split()
                    data.append(split[self.skipCols:])
                    self.skippedCols.append(split[:self.skipCols])

                if self.skipCols != 0:
                    self.skippedCols = zip(*self.skippedCols)

                data = numpy.array(data, dtype=numpy.float)

#                data = numpy.array(self.rx.findall(self.file.read()), numpy.float)
                data = data if self.transpose else data.transpose()


            if time.time() - t0 > 10.0:
               self.Error("TIMEOUT: File was empty for too long...")
               self.file.close()
               self.Exit()
               return
               
        self.file.close()
        
        return dataGenerator(data)
    
    def unpackNumpyBin(self, binFile):
        
        binFile.seek(0)
        try:
            data = numpy.load(binFile)
            return data
        except:
            return []
    
    def unpackBinFile(self, binFile):
 
        data = numpy.fromfile(binFile, dtype=numpy.float64)
  
        m = self.Ncols
        n = int(data.size/m)
        
        data.resize(n, m)
        
        return data if not self.transpose else data.transpose()
            
    
    def unpackArmaMatBin(self, armaFile):   

        armaFormat =  armaFile.readline()

        dims = tuple([int(d) for d in armaFile.readline().strip().split()])
      
        if "IS004" in armaFormat:
            dtype = numpy.int32
        else:
            dtype = numpy.float64
        
        
        
        data = numpy.fromfile(armaFile, dtype=dtype).transpose()

        if len(dims) == 2:

            _data = numpy.zeros(dims)            
            for i in range(dims[1]):
                _data[:, i] = data[i*dims[0]:(i+1)*dims[0]]
                
            data = _data

        if self.transpose:
            data.resize(dims)
            data = data.transpose()
        else:
            data.resize(dims)
        
        return data
    
    def set_figures(self):
        
        if self.stack not in ["H", "V"]:
            self.Error("Invalid stack argument %s. (Choose either H or V)" % self.stack)
            self.Exit()
        
        s = ""
        i = 0
        self.figures = []
        
        if self.makeGif:
            self.savedImages = []
        
        if not self.figMap:
            self.figMap = {"figure": ["subfigure"]}
        
        for fig in self.figMap.keys():
            s += "self.%s = plab.figure(); " % fig
            s += "self.i%s = self.add_figure(self.%s); " % (fig, fig)
        
            subFigs = self.figMap[fig]
            
            if self.makeGif:
                self.savedImages.append([])

            if type(subFigs) is str:
                subFigs = [subFigs]
            
            nFigs = len(subFigs)
            
            for j in range(nFigs):
                
                
                if self.stack == "V":
                    s += "self.%s = self.%s.add_subplot(%d, 1, %d); " % (subFigs[j], fig, nFigs, j+1)
                else:
                    s += "self.%s = self.%s.add_subplot(1, %d, %d); " % (subFigs[j], fig, nFigs, j+1)
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
                    self.Exit()
                    return          


    def Error(self, s):
        if self.useGUI:
            self.parent.parent.terminalTracker("DCViz", s)
        else:
            print "[%s] %s" % ("DCViz".center(10), s)


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
            if not self.plotted:
                self.show(drawOnly = self.toFile)
            else:
                self.show(drawOnly=True)
        else:
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

    def initFamily(self):

        self.originalFilename = os.path.split(self.filepath)[-1]
        
        if self.isFamilyMember and self.loadSequential and self.smartIncrement and self.getNumberForSort:
        
            familyMembers = self.get_family()            
            
            if self.getNumberForSort:
                familyMembers = sorted(familyMembers, key=self.getNumberForSort)
                
                try:
                    inc = self.getNumberForSort(familyMembers[1]) - self.getNumberForSort(familyMembers[0])    
                except IndexError:
                    print "Only one familiyMember present in directory.. set 'smartIncrement = false' will remove this error"
                    inc = 1
            self.ziggyMagicNumber /= inc
            
            if self.ziggyMagicNumber < 1:
                self.ziggyMagicNumber = 1
                
    def mainloop(self):
        
        if not self.armaBin:
            breakMe = self.waitForGreenLight()
            
            if breakMe:
                return

        if self.hugifyFonts:
            self.hugify()

        self.manageFigures()
        self.initFamily()

        while self.shouldReplot():

            self.clear()
           
            data = self.get_data(setUpFamily=self.isFamilyMember)

            self.plot(data)  
            self.showFigures()
            self.plotted = True

            if self.dynamic:
                
                if self.makeGif:
                    self.prepGif()
                else:
                    self.sleep() 
                    
            else:
                if not self.useGUI and not self.toFile and not self.makeGif:
                    self.stall()
                break
                
        if not self.useGUI:
            if self.toFile:
                self.saveFigs()
            elif self.makeGif:
                self.makeGifs()
            self.close()

    def prepGif(self):
        self.saveFigs()
    
    def makeGifs(self):
        
        
            
        path, fname = os.path.split(self.filepath)
        
        dirname = "DCViz_out"
        dirpath = pjoin(path, dirname)
        
        thisDir = os.getcwd()
        os.chdir(dirpath)

        for i in range(len(zip(*self.figures)[0])):
            
            if len(self.savedImages[i]) < 2:
                print "Making a Gif requires more than 1 image."
                print "images: ", self.savedImages
                continue

            imgs = " ".join(self.savedImages[i])
            gifName = ".".join(fname.split(".")[0:-1]) + "_" + str(i) + ".gif"

            print "Making gif %s (%d/%d) out of %s-%s" % (gifName, i+1, len(zip(*self.figures)[0]), self.savedImages[i][0], self.savedImages[i][-1])
            os.system("convert -delay %g %s -loop %g %s" % (self.delay, imgs, self.gifLoopDelay, gifName))            
        
        os.chdir(thisDir)
                   
           
    def stall(self):
        raw_input("[%s] Press any key to exit" % "DCViz".center(10))

    def rawDataAPI(self, data):
        if self.dynamic:
            print "[%s] Dynamic mode not supported for direct push mode." % "DCViz".center(10)
        
        if self.filepath is None:
            if self.toFile is True:
                print "[%s] A filename needs to be supplied in order to save figures to file." % "DCViz".center(10)
        else:
            self.filepath = self.filepath.strip(".png")
        data = dataGenerator(data)
        
        self.manageFigures()
        self.plot(data)
        self.showFigures()
        
        if self.toFile:
            self.saveFigs()
        else:
            self.stall()
            
        self.close()
            
        
    def load_sample(self):
        
        self.reload()
  
        sample = self.file.read()
        
        if not sample:
            self.Error("No data in file...")
            return "red"
            
        if self.fileBin:
            return "green"
            
        self.skipRows, self.Ncols = self.sniffer(sample) 
        
        if (self.skipRows == 0):
            self.skipRows = None
        
        self.file.seek(0)        

        if self.skipRows:
            
            self.skippedRows = []
            for i in range(self.skipRows):
                self.skippedRows.append(self.file.readline().strip())
      
      
        return "green"
        
    
    def sniffer(self, sample):
        
        sampleList = [row.split() for row in sample.split("\n")]
 
        #If user has specified the number of rows to skip
        if self.skipRows is not None:
            return self.skipRows, len(sampleList[0])

        nLast = len(sampleList[-1])
        i = len(sampleList) - 1

        if nLast != len(sampleList[-2]):
            nLast = len(sampleList[-2])
            i -= 1
        
        while len(sampleList[i]) == nLast and i >= 0:
            i -= 1

        #(nRows to skip = i+1, nCols = nLast)
        return i+1, nLast    
    
    def recursiveSafeLoad(self, depth, maxDepth):
        
        N = 0        
        while not os.path.exists(self.filepath) and N < 10:
                time.sleep(0.1)
                N += 1
        
        try:
            self.file = open(self.filepath , "r")
        except:
            if (depth == maxDepth):
                self.Error("Unable to load file.")
                sys.exit(1) 
                
            self.recursiveSafeLoad(depth + 1, maxDepth)
        
    
    def reload(self):

        self.path, self.filename = os.path.split(self.filepath)
        
        if self.file:
            if not self.file.closed: 
                self.file.close()

        self.recursiveSafeLoad(0, 3)
   
        self.skippedRows = []
        if not self.armaBin and self.skipRows is not None:        
            for i in range(self.skipRows):
                self.skippedRows.append(self.file.readline().strip())

                        
        if self.fileBin and self.binaryHeaderBitSizes:
            self.binaryHeader = []
            offset = 0
            for i, size in enumerate(self.binaryHeaderBitSizes):
                
                if self.binaryHeaderTypes:
                    bintype = self.binaryHeaderTypes[i]
                else:
                    bintype = self.defaultBinaryHeaderTypes[size]

                try:
                    self.binaryHeader.append(struct.unpack(bintype, self.file.read(size))[0])
                except struct.error:
                    return

                offset += size
            
            if self.nColsFromHeaderLoc:
                self.Ncols = self.binaryHeader[self.nColsFromHeaderLoc]


    def saveFigs(self):
            
        path, fname = os.path.split(self.filepath)
        
        dirname = "DCViz_out"
        dirpath = pjoin(path, dirname)
        
        if not os.path.exists(dirpath):
            os.mkdir(dirpath)
        
        i = 0
        for fig in zip(*self.figures)[0]:
            
            figname = ".".join(fname.split(".")[0:-1]) + "_" + str(i) + ".png"
            figpath = pjoin(dirpath, figname)
            fig.savefig(figpath, transparent=self.transparent)

            if self.makeGif:
                self.savedImages[i].append(figname)

            i += 1
        
        print "[%s] Figure(s) successfully saved." % "DCViz".center(10)
        
    def hugify(self):
        rcParams['font.size'] = self.fontSize
        rcParams['xtick.labelsize'] = self.labelSize
        rcParams['ytick.labelsize'] = self.labelSize
        
#        rcParams['ytick.major.size'] = self.tickSize
#        rcParams['ytick.major.width'] = self.tickSize
#        
#        rcParams['xtick.major.size'] = self.tickSize
#        rcParams['xtick.major.width'] = self.tickSize
        
        rcParams['axes.labelsize'] = self.labelSize
        rcParams['axes.titlesize'] = self.labelSize

    def add_figure(self, fig):
        self.figures.append([fig])
        self.nFig += 1
        return len(self.figures) - 1
        
    def add_subfigure(self, subfig, i):
        self.figures[i].append(subfig)
        return len(self.figures[i]) - 1
        
    def show(self, drawOnly=False):
        for fig in self.figures:
#            try:
            fig[0].canvas.draw()
#            except:
#                raise OSError("Unable to draw canvas! Missing dvips drivers? (sudo apt-get install dvips-fontdata-n2bk)\nAlso check all your labels etc. for weird symbols or latex-expressions without assigned $-signs.")
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
    
    def Exit(self):
        
        if self.useGUI:
            self.parent.closeEvent(None)
        else:
            sys.exit()
        
    
    def plot(self, data):
        raise NotImplementedError("Plot function must be implemented for a DCViz instance.")
    
    
