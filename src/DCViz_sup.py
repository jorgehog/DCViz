# -*- coding: utf-8 -*-

import re, numpy, time, signal, os, sys, struct, itertools

import matplotlib.pylab as plab

from matplotlib import rcParams
from os.path import join as pjoin
from random import shuffle

class DCVizLoader(object):

    def __init__(self):
        self.plotter = None
        self.skippedCols = []
        self.skippedRows = []

    def get_metadata(self):
        return None

    def set_plotter_instance(self, instance):
        self.plotter = instance

    def load(self, file):
        raise RuntimeError("Load method is not implemented")

    def load_from_path(self, path):

        if not os.path.exists(path):
            raise RuntimeError("Invalid path: %s" % path)

        with open(path, 'r') as f:
            return self.load(f)

class RawAscii(DCVizLoader):

    def __init__(self, dtype=numpy.float):
        self.dtype = dtype

        super(RawAscii, self).__init__()

    def sniffer(self, sample):

        sampleList = [row.split() for row in sample.split("\n")]

        #If user has specified the number of rows to skip
        if self.plotter.skipRows is not None:
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

    def get_dims(self, file):
        self.plotter.skipRows, self.plotter.Ncols = self.sniffer(self.plotter.sample)

        if (self.plotter.skipRows == 0):
            self.plotter.skipRows = None


    def load(self, file):
        self.get_dims(file)

        data = []

        for i in range(self.plotter.skipRows):
            self.skippedRows.append(file.readline())

        for line in file:
            split = line.split()
            data.append(split[self.plotter.skipCols:])
            self.skippedCols.append(split[:self.plotter.skipCols])

        if self.plotter.skipCols != 0:
            self.skippedCols = zip(*self.skippedCols)

        return numpy.array(data, dtype=self.dtype)

    def get_metadata(self):
        return self.skippedCols, self.skippedRows

class Numpy(DCVizLoader):
    def load(self, file):
        return numpy.load(file)

class Armadillo(DCVizLoader):

    def load(self, file):
        armaFormat = file.readline()

        dims = tuple([int(d) for d in file.readline().strip().split()])

        if 0 in dims:
            print "Zero dimension array loaded.", dims
            return []

        if "IS004" in armaFormat:
            dtype = numpy.int32
        else:
            dtype = numpy.float64

        data = numpy.fromfile(file, dtype=dtype).transpose()

        if len(dims) == 2:

            _data = numpy.zeros(dims)
            for i in range(dims[1]):
                _data[:, i] = data[i*dims[0]:(i+1)*dims[0]]

            data = _data

        data.resize(dims)

        return data

class BinaryWithHeader(DCVizLoader):

    default_binary_header_types = {4: 'i', 8: 'd'}

    def __init__(self, binary_header_bit_sizes, n_cols_header_index=None, binary_header_types=None):

        self.binary_header_bit_sizes = binary_header_bit_sizes
        self.binary_header_types = binary_header_types
        self.n_cols_header_index=n_cols_header_index

        self.binary_header = []

        super(BinaryWithHeader, self).__init__()

    def get_metadata(self):
        return self.binary_header

    def load_header(self, file):

        self.binary_header = []

        offset = 0
        for i, size in enumerate(self.binary_header_bit_sizes):

            if self.binary_header_types:
                bintype = self.binary_header_types[i]
            else:
                bintype = self.default_binary_header_types[size]

            self.binary_header.append(struct.unpack(bintype, file.read(size))[0])

            offset += size


    def load(self, file):

        self.load_header(file)

        data = numpy.fromfile(file, dtype=numpy.float64)

        if self.n_cols_header_index:
            m = self.binary_header[self.n_cols_header_index]
        else:
            m = self.plotter.Ncols

        n = int(data.size/m)

        data.resize(n, m)

        return data

class Binary(BinaryWithHeader):

    def __init__(self):
        super(Binary, self).__init__([])

class Ignis(BinaryWithHeader):

    def __init__(self):
        super(Ignis, self).__init__([4,4],n_cols_header_index=1)



class DCVizPlotter:

    subfigure = None
    figMap = {}
    
    nametag = None

    loader = None
    loader_ending_map = {"arma" : Armadillo,
                         "npy" : Numpy,
                         "ign" : Ignis}

    path = None
    filename = None
    filepath = None

    Ncols = None
    sample = None
    transpose = False
    
    isFamilyMember = False
    familyName = None
    familyFileNames = []
    familyHome = None
    familyHead = None
    loadLatest = False
    loadSequential = False
    getNumberForSort = lambda _s, x : int(re.findall(_s.nametag, x)[0])
    ziggyMagicNumber = 1
    smartIncrement = True
    originalFilename = None

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
    labelSize= 25
    ticklabelSize=20
    fontSize = 20  # Only invoked by hugifyFonts = True
    tickSize = 2
    
    fig_size = None
    relative_fig_size = None

    tight = True

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
        self.argv = []
        
        self.filepath = filepath
        self.file = None
        
        self.toFile = toFile

        if self.loadSequential:
            self.nextInLine = 0
        
        if not (threaded or useGUI) and dynamic:
            print "[  DCViz  ]", "Interrupt dynamic mode with CTRL+C"
            signal.signal(signal.SIGINT, self.signal_handler)

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

    def get_family_index_from_name(self, name):
        return self.familyFileNames.index(name)

    def get_family_member_data(self, data, unique_name):
        pattern = re.sub("\(.*\)", unique_name, self.nametag)

        for family_name in self.familyFileNames:

            if re.findall(pattern, family_name):
                return data[self.get_family_index_from_name(family_name)]

        raise RuntimeError("Name '%s' not found in family." % unique_name)

    def get_family(self):
        
        familyHome, self.familyHead = os.path.split(self.filepath)
            
        if not familyHome:

            if self.useGUI:
                self.Error("Single file name given as path. DCViz needs the absolute file path to work with families.")
                self.Exit()
            else:
                familyHome = os.getcwd()

        familyNames = [name for name in os.listdir(familyHome)\
                        if re.findall(self.nametag, name) and os.path.exists(pjoin(familyHome, name)) and "tmp" not in name]
       
        familyMembers = sorted([pjoin(familyHome, name) for name in familyNames])     

        self.familyHome = familyHome

        return familyMembers

    def dictify(self, data):

        _data = {}        
        for fdata, name in zip(data, self.familyFileNames):
            _data[name] = fdata
            
        data = _data

    def get_nametag_match(self, name):

        match = re.findall(self.nametag, name)

        if not match:
            return ""

        return match[0]

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

                self.family_loader_data = []

                data = [0]*N
                self.familyFileNames = [0]*len(data)     
                
                for i in range(N):
                 
                    self.file = None
                    self.filepath = familyMembers[i]
                
                    self.familyFileNames[i] = os.path.basename(familyMembers[i])
                    
                    data[i] = self.get_data(setUpFamily=False)
                    
                    self.family_loader_data.append(self.loader.get_metadata())

                    
            self.file.close()
            return data

        if self.loader is None:

            #If we have a name like .arma og .npy or .ign then we automatically select the loader
            ending = self.filepath.split(".")[-1]

            if ending in self.loader_ending_map.keys():
                self.loader = self.loader_ending_map[ending]()
            else:
                self.loader = RawAscii()

        data = []

        #attempt to reload file untill data is found.
        t0 = time.time()
        while len(data) == 0:

            self.reload()



            # if self.armaBin:
            #     loader = Armadillo()
            #     # data = self.unpackArmaMatBin(self.file)
            # elif self.fileBin:
            #     loader = BinaryWithHeader(self.binaryHeaderBitSizes, self.nColsFromHeaderLoc, self.binaryHeaderTypes)
            #     # data = self.unpackBinFile(self.file)
            # elif self.numpyBin:
            #     loader = Numpy()
            #     # data = self.unpackNumpyBin(self.file)
            # else:
            #     loader = RawAscii()
            #     # data = []
            #     #
            #     # for line in self.file:
            #     #     split = line.split()
            #     #     data.append(split[self.skipCols:])
            #     #     self.skippedCols.append(split[:self.skipCols])
            #     #
            #     # if self.skipCols != 0:
            #     #     self.skippedCols = zip(*self.skippedCols)
            #     #
            #     # data = numpy.array(data, dtype=numpy.float)

            self.loader.set_plotter_instance(self)
            data = self.loader.load(self.file)
            data = data if not self.transpose else data.transpose()


            if time.time() - t0 > 10.0:
               self.Error("TIMEOUT: File was empty for too long...")
               self.file.close()
               self.Exit()
               return
               
        self.file.close()

        #In case of vectors we get rid of the wrapping array
        if data.size != 1 and len(data.shape) != 1:
            if len(data) == 1:
                data = data[0]
            elif data.shape[1] == 1:
                data = data[:, 0]

        return data
    
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
        
        return data
            
    
    def unpackArmaMatBin(self, armaFile):   

        armaFormat =  armaFile.readline()

        dims = tuple([int(d) for d in armaFile.readline().strip().split()])

        if 0 in dims:
            print "Zero dimension array loaded.", dims
            return []

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

        data.resize(dims)

        return data
    
    def set_figures(self):
        
        if self.stack not in ["H", "V"]:
            self.Error("Invalid stack argument %s. (Choose either H or V)" % self.stack)
            self.Exit()
        
        s = ""
        i = 0
        self.figures = []
        self.figure_names = []
        
        if self.makeGif:
            self.savedImages = []
        
        if not self.figMap:
            self.figMap = {"figure": ["subfigure"]}
        
        for fig in self.figMap.keys():
            s += "self.%s = plab.figure(); " % fig
            s += "self.i%s = self.add_figure(self.%s, '%s'); " % (fig, fig, fig)
        
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
                
    def mainloop(self, argv=[]):
        self.argv = argv

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
  
        self.sample = self.file.read()
        
        if not self.sample:
            self.Error("No data in file...")
            return "red"
        else:
            return "green"

        self.file.seek(0)

    #
    # def sniffer(self, sample):
    #
    #     sampleList = [row.split() for row in sample.split("\n")]
    #
    #     #If user has specified the number of rows to skip
    #     if self.skipRows is not None:
    #         return self.skipRows, len(sampleList[0])
    #
    #     nLast = len(sampleList[-1])
    #     i = len(sampleList) - 1
    #
    #     if nLast != len(sampleList[-2]):
    #         nLast = len(sampleList[-2])
    #         i -= 1
    #
    #     while len(sampleList[i]) == nLast and i >= 0:
    #         i -= 1
    #
    #     #(nRows to skip = i+1, nCols = nLast)
    #     return i+1, nLast
    
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

        # if self.fileBin and self.binaryHeaderBitSizes:
        #     self.binaryHeader = []
        #     offset = 0
        #     for i, size in enumerate(self.binaryHeaderBitSizes):
        #
        #         if self.binaryHeaderTypes:
        #             bintype = self.binaryHeaderTypes[i]
        #         else:
        #             bintype = self.defaultBinaryHeaderTypes[size]
        #
        #         try:
        #             self.binaryHeader.append(struct.unpack(bintype, self.file.read(size))[0])
        #         except struct.error:
        #             return
        #
        #         offset += size
        #
        #     if self.nColsFromHeaderLoc:
        #         self.Ncols = self.binaryHeader[self.nColsFromHeaderLoc]


    def saveFigs(self):

        if self.isFamilyMember:
            filepath = pjoin(self.familyHome, self.familyHead)
        else:
            filepath = self.filepath

        path, fname = os.path.split(filepath)
        
        dirname = "DCViz_out"
        dirpath = pjoin(path, dirname)
        
        if not os.path.exists(dirpath):
            os.mkdir(dirpath)
        
        i = 0
        for fig in zip(*self.figures)[0]:
            
            figname = ".".join(fname.split(".")[0:-1]) + "_" + self.figure_names[i] + ".png"
            figpath = pjoin(dirpath, figname)
            fig.savefig(figpath, transparent=self.transparent)

            if self.makeGif:
                self.savedImages[i].append(figname)

            i += 1
        
        print "[%s] Figure(s) successfully saved." % "DCViz".center(10)
        
    def hugify(self):
        rcParams['font.size'] = self.fontSize
        rcParams['xtick.labelsize'] = self.ticklabelSize
        rcParams['ytick.labelsize'] = self.ticklabelSize
        
#        rcParams['ytick.major.size'] = self.tickSize
#        rcParams['ytick.major.width'] = self.tickSize
#        
#        rcParams['xtick.major.size'] = self.tickSize
#        rcParams['xtick.major.width'] = self.tickSize
        
        rcParams['axes.labelsize'] = self.labelSize
        rcParams['axes.titlesize'] = self.labelSize

    def add_figure(self, fig, figname):
        self.figures.append([fig])
        self.figure_names.append(figname)
        self.nFig += 1
        return len(self.figures) - 1
        
    def add_subfigure(self, subfig, i):
        self.figures[i].append(subfig)
        return len(self.figures[i]) - 1
        
    def show(self, drawOnly=False):
        for fig in self.figures:

            if self.tight:
                fig[0].tight_layout()

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
