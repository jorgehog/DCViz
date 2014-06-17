
# -*- coding: utf-8 -*-

import sys, re, os, inspect, datetime
from re import findall as find

try:
    import numpy
    import numpy as np
    from numpy import *
except:
    print "\n"
    print "You need numpy in order to run this library!"
    print "sudo apt-get install python-numpy-dev"
    print "\n" 
    sys.exit(1)


#==============================================================================
# Setting font-style to latex if possible
#==============================================================================

from matplotlib import rc
from matplotlib import rcParams

try:
    rc('text', usetex=True)
    rc('font', family='serif')
    
    try:
        from matplotlib import pylab
        from matplotlib.pylab import *
    except:
        print "\nYou need matplotlib to use DCViz."
        print "sudo apt-get install python-matplotlib"
        print "\n"
        sys.exit(1)
    
    pylab.title("$\LaTeX$")
    pylab.draw()
    pylab.clf()
    
except:
    print "Neccessary latex packages not installed. Disabling latex support."
    print "Got Latex? Make sure you have the DVIPS package from"
    print "\nsudo apt-get install dvips-fontdata-n2bk\n\n"

    rc('text', usetex=False)   


#==============================================================================
# Additional MPL imports goes here
#==============================================================================

from matplotlib import pylab, colors, ticker, cm
from mpl_toolkits.mplot3d import Axes3D


#==============================================================================
# Including DCViz superclass and additional tools
#==============================================================================
classes_thisDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(classes_thisDir)

from DCViz_sup import DCVizPlotter, dataGenerator


#==============================================================================
# User-defined classes:
#==============================================================================


class myTestClass(DCVizPlotter):
    nametag =  'testcase\d\.dat' #filename with regex support
    
    #1 figure with 1 subfigure
    figMap = {'fig1': ['subfig1']}
    
    #skip first row. (the function __str__ is printed here)
    skipRows = 1    
    
    def plot(self, data):
        column1 = data[0]

        self.subfig1.set_title('I have $\LaTeX$ support!')
              
        self.subfig1.set_ylim([-1,1])
        self.subfig1.set_xlim([0, 50])
          
        self.subfig1.plot(column1)


class results(DCVizPlotter):
    nametag = "msd\.dat|vacf\.dat"
    isFamilyMember = True
    
    def plot(self, data):
        
        for i, fdata in enumerate(data):
            t, vr, d = fdata
            self.subfigure.plot(t, d, label=self.familyFileNames[i].strip(".dat"))
        
        pstar = 0.2
        Tstar = 0.5
        self.subfigure.plot(t, numpy.ones(len(t)) + 10**(0.05 + 0.07*pstar - (1.04 + 0.1*pstar)/Tstar))
        self.subfigure.axes.set_ylabel("D")
        self.subfigure.legend()

class myTestClassFamily(DCVizPlotter):
    nametag =  'testcaseFamily\d\.dat' #filename with regex support
    
    #1 figure with 3 subfigures
    figMap = {'fig1': ['subfig1', 'subfig2', 'subfig3']}
    
    #skip first row. (the function __str__ is printed here)
    skipRows = 1    

    #Using this flag will read all the files matching the nametag
    #(in the same folder.) and make them aviable in the data arg    
    isFamilyMember = True
    familyName = "testcase"
    
    def plot(self, data):
        
        #figures[0] is 'fig1' figures. the 0'th element is the
        #self.fig1 itself. Subfigures are always index [1:]
        mainFig = self.figures[0][0]  
        mainFig.suptitle('I have $\LaTeX$ support!')        
        subfigs = self.figures[0][1:]
    
        #Notice we plot fileData.data and not fileData alone.
        #The dataGenerator class is used to speed up file reading;
        #looping over family members and directly plotting means
        #we send a dataGenerator instance to matplotlib.
        #in order to get the numpy object, we send the data.
        #Alternatively, we could send data[:]
        for subfig, fileData in zip(subfigs, data):
            subfig.plot(fileData.data)
            subfig.set_ylim([-1,1])
        

class standardBinaryArmaVec(DCVizPlotter):
    
    nametag = "binaryArmaVec\.arma$"

    armaBin = True
    
    figMap = {'fig' : 'subfigure', 'fig2': 'sf'}
    def plot(self, data):

        N = len(data)        
        r = numpy.linspace(0.1, 6, N)
        self.subfigure.plot(r, numpy.cumsum(data.data)*(6-0.1)/(N-1), 'b-^');
        self.sf.plot(r, -data.data, 'b-^')
        self.sf.axes.set_xlabel("R")
        self.subfigure.axes.set_xlabel("R")
        self.sf.axes.set_ylabel("F")
        self.subfigure.axes.set_ylabel("-cumsum(F)")
        
from scipy.stats import linregress

class forces1D(DCVizPlotter):
    
    nametag = "forces\.arma"
    
    armaBin = True
    
    figMap = {"fig" : ("subfigure", "figure")}
    def plot(self, data):
        
        all_dr = linspace(0.9, 1, len(data));      
        
        dr = (1 - all_dr)[::-1]

        data.data = data.data[::-1]
        
        self.subfigure.plot(dr, data.data)
        self.subfigure.set_xlabel("compression")
        self.subfigure.set_ylabel("wall force")
        
        self.figure.plot(dr, log(data.data))
        self.figure.set_xlabel("compression")
        self.figure.set_ylabel("log(wall force)")
        
        b = len(data.data)/20
        slope, intercept, rV, pV, stdErr = linregress(dr[b:], log(data.data)[b:][:, 0])
        
        print dr[b:]
        self.figure.plot(dr, slope*dr + intercept, "k", label="%s" % slope)
        pylab.legend(loc=0)
        
        
class testStuff(DCVizPlotter):
    
    nametag = "testStuff\.arma"
    
    armaBin = True
    
    def plot(self, data):
        
        x, y, = data
        
        print x.shape, y.shape
        self.subfigure.scatter(x, y, s = 1000)
      
class KMC_1D(DCVizPlotter):

    nametag = "1DKMC\.npy"
    
    numpyBin = True
    
    figMap = {"fig": ["zFig", "Efig", "Rfig", "Afig"]}
    
#    relative_fig_size = [2, 2]    
    
    def plot(self, data):

        z, E, Rl, Rr, A1, A2 = data
        
        A = numpy.concatenate((A1, A2))        
        
        LX, LY = pylab.rcParams['figure.figsize']
        LY /= 5
                
        

        self.zFig.bar(numpy.arange(len(z)), z, linewidth = 0, width=1)
#        self.zFig.set_ylim(0, LY/LX*len(z))
        self.Efig.plot(E, "b*")
        self.Efig.set_ybound(0)
        self.Rfig.plot(Rl, "b^")
        self.Rfig.plot(Rr, "r*")
        self.Rfig.set_ybound(0)
        self.Rfig.set_xlim([-1, len(z)])
        self.Afig.plot(A, "bo")
        self.Afig.set_ybound(0)
      
class concentrations(DCVizPlotter):
    
    nametag = 'concOut.+?\.arma'
    
#    figMap = {'fig' : ['sf_c', 'sf_CO2', 'sf_Ca', 'sf_mass', 'sf_cl', 'sf_ph']}
    figMap = {'fig' : ['sf_c'], 'fig2' : 'sf_CO2', 'fig3' : 'sf_Ca', 'fig6' : 'sf_mass', 'fig7' : 'sf_cl', 'fig8' : ['sf_ph']}
    species = ['c', 'CO2', 'Ca', 'mass', 'cl', 'ph']
    
    ylim = {'c' : [0], 
            'CO2' : [0], 
            'Ca' : [], 
            'mass': [0], 
            'cl': [0], 
            'ph': [1, 14]}

    ylab = {'c' : [], 
            'CO2' : [], 
            'Ca' : [], 
            'mass': 'Density CaCO3 [mg/cm$^2$]', 
            'cl': [], 
            'ph': 'pH'}
            
    nightmareAxes = ['c', 'CO2', 'Ca']

    armaBin = True
    transposed = True    
    
    isFamilyMember = True
    loadLatest = True
#    loadSequential = True
    
    getNumberForSort = lambda _s, a: int(find("concOut(\d+?)\.arma", a)[0])
    
    stack = "H"
    
    ziggyMagicNumber = 10000
    
    gifLoopDelay = 5
    
    txt = None
    lsize = 20
    tsize = 20
    
    def plot(self, data):
        
        dt = 0.01       
        
        for i, spec in enumerate(self.species):
            subfigure = eval('self.sf_' + spec)
            
            subfigure.plot(data[i], 'b-*', markersize=20, markeredgewidth=1, markeredgecolor='k')
            
            subfigure.set_xlabel('Cell \#', size=self.lsize)
            
            if spec in self.nightmareAxes:
                formatter = ticker.ScalarFormatter(useOffset=True)
                formatter.set_scientific(True)
                formatter.set_powerlimits([-3, 2])
                subfigure.axes.get_yaxis().set_major_formatter(formatter)
                
                subfigure.axes.get_yaxis().offsetText.set_size(self.tsize)

                
            for tick in subfigure.axes.yaxis.get_major_ticks():
                tick.label.set_fontsize(self.tsize)  
                
            subfigure.axes.get_xaxis().set_major_locator(pylab.MaxNLocator(integer=True));
            
            for tick in subfigure.axes.xaxis.get_major_ticks():
                tick.label.set_fontsize(self.tsize)             
                               
                             
            
            if self.ylim[spec]:
                if len(self.ylim[spec]) == 1:
                    subfigure.axes.set_ybound(self.ylim[spec][0])
                else:
                    subfigure.axes.set_ylim(self.ylim[spec][0])
                    
            if self.ylab[spec]:
                subfigure.set_title(self.ylab[spec], size=self.tsize)
            else:
                subfigure.set_title("[%s]" % spec, size=self.tsize)                
                
  
#            
        __N = self.getNumberForSort(self.filepath)
        s = str(datetime.timedelta(seconds=__N*dt)).split(".")[0]
#        s = str(__N)
        print s
#        if not self.txt:        
#            self.txt = self.fig.suptitle(s)
#        else:
#            self.txt.set_text(s)
##        
#        self.sf_c.set_title("[c]")
#        self.sf_c.set_ylabel("conc. [M/l]")
                
        

class phreeqc_si(DCVizPlotter):
    
    nametag = 'phreeqc.*\.dat'

    skipRows = 1
    skipCols = 5
    
#    figMap = {"fig": ['subfigure', 'molfig',  'hfig', 'phfig']}
    figMap = {"fig": ['subfigure', 'phfig']}
    
    nTotals = 4
    
    stack = "V"
    
    hugifyFonts = True    
    fontSize=15    
    labelSize=40
    tickSize=1
    
    coPrIt = 1.0    
    
    def transName(self, name):
        return name.replace("Brucite", "Mg(OH)2")
    
    def plot(self, data):
        
        
        lol = self.getRandomStyles()
    
        lol = ['r-*', 'g-+', 'b-^', 'b--x', 'r-o', 'k-x', 'g--^']
        step, pH, pe, mg, F2, F3, H = data.data[:(3 + self.nTotals)]
        step[0] = 0
        step = numpy.array(step)
        step *= self.coPrIt

        self.phfig.plot(step, pH, "k-x", label='pH')
#        self.phfig.plot(step, -numpy.log10(F2), label='pF++')
#        self.phfig.plot(step, -numpy.log10(F3), label='pF+++')
#        self.phfig.plot(step, -numpy.log10(H), label='pH2')
        self.phfig.plot(step, -numpy.log10(mg), 'r-x', label='pMg')
        self.phfig.set_xlabel("Added CO2 [mol]")
        self.phfig.set_ylabel("pX")       
        self.phfig.ticklabel_format(useOffset=False, axis='y')
        self.phfig.legend()
        self.phfig.set_xbound(1)
        
        k = self.nTotals        
        
#        for k, _m in enumerate(data.data[3:3+self.nTotals]):
#            name = self.skippedRows[0].split()[8 + k]
#            
#            print _m, name
#            if "H(0)" not in name:
#                self.molfig.plot(step, _m, lol[3*k], label=name)
#            else:
#                self.hfig.plot(step, _m)
#                self.hfig.set_ylabel('Conc H2')
#            
#        k += 1
        j = 0
        for i, phase in enumerate(data.data[(4+k):][::2]):

            name = self.skippedRows[0].split()[8 + k + 2*i]    

            if len(numpy.where(phase != 0)[0]) == 0:        
                continue
     
            self.subfigure.plot(step, phase, lol[j], label=self.transName(name))
            j+=1

#        self.subfigure.set_xlabel("Reaction step")

#        self.molfig.legend()
#        self.molfig.set_ylabel("Total concentration")        
#        
        self.subfigure.set_ylabel("Change [mol]")
        self.subfigure.legend(loc=0)
        self.subfigure.set_xbound(1)

class molecules(DCVizPlotter):
    
    nametag = "dist\_out\_Molecule.+?arma3D"
    
    armaBin = True    

    
    def plot(self, data):
        
        data = numpy.sum(data.data, axis=2)
        
        n = data.shape[0]

        data = data[n/4:3.2*n/4, n/3:3*n/4]        
        
        im = self.subfigure.imshow(data, cmap=pylab.cm.spectral, interpolation="lanczos");
        self.figure.colorbar(im)


#==============================================================================
# avconv -f image2 -y -r 24 -i mdPos%05d_0.png -b:v 1000k simple2D_LJ.avi
#==============================================================================

class forPress(DCVizPlotter):
    
    nametag = "FOR_PRESS\.tmp"
    
    figMap = {"fig": ["cl", "qm"]}
    
    cl_v = [0]
    qm_v = []    
    
    def superAppend(self, *args):
        
        for a in args:
            self.qm_v.append(a)
    
    def plot(self, data):
    
        self.qm.axes.set_ylabel("P(M)")        
        self.cl.axes.set_ylabel("P(M)")
        self.qm.axes.set_xlabel("Measurement (M)")
        
        self.cl_v.append(self.cl_v[0])
        self.superAppend(*list(numpy.random.normal(size=len(self.cl_v))))
        
        if len(self.cl_v) == 2:
            return
        
        self.cl.hist(self.cl_v, 20, normed=True)
        self.qm.hist(self.qm_v, 20, normed=True)

        self.cl.axes.set_xlim([-3, 3])
        self.qm.axes.set_xlim([-3, 3])
        
        self.qm.axes.set_ylim([0, 1])        
        self.cl.axes.set_ylim([0, 1])
        
        path, name = os.path.split(self.filepath)
        _file = os.path.join(path, "DCViz_out", "FOR_PRESS_0.png")
        if os.path.exists(_file):
            os.rename(_file, _file.replace("_0", "_" + str(len(self.cl_v)-2).rjust(5, "0")))



class virials(DCVizPlotter):
    
    nametag = "w\_t\_HO\_COL\_N(\d+)\.dat"
    
    isFamilyMember = True
    
#    figMap = {"fig" : ("subfigure", "phasefigure")}
    
    transpose = True    
    
    def plot(self, data):
        
        for i, _data in enumerate(data):
            N = re.findall(self.nametag, self.familyFileNames[i])[0]
            

            w, alpha, beta, E, T, vho, vcol, r, r2, rij, err_E, err_T, err_vho, err_vcol, err_r, err_r2, err_rij = _data
        
            if vcol.sum() != vcol.sum():
                vcol = 0
    
            V = vho + vcol
            
            sigmaGAMMA = (err_T*V + T*(err_vho + err_vcol))/(err_T**2 - T**2)
            
            I = numpy.where(V > 1E-10)
            

            l = "N=%s" % N              
            
            print l
            for o, E, _err_E in zip(w, T+V, err_E):
                print "%g   %g    %g"  % (o, E, _err_E)
            print "------------------------------"
            
            GAMMA = T[I]/V[I]
            
#            DGAMMA = numpy.diff(GAMMA)            
            
            
            self.subfigure.plot(w[I], GAMMA, label=l, marker='^', linestyle="--")
#            self.phasefigure.plot(w[I][:-1], DGAMMA, label=l, marker='^', linestyle="--")            
            
            
#        self.phasefigure.plot(T[I], T[I], "k--", label="NoCoulomb")

        self.subfigure.legend(loc=4)
        self.subfigure.set_xlabel("$\omega$")
        self.subfigure.set_ylabel("T/V")

#        self.phasefigure.legend()        
#        self.phasefigure.set_ylabel("D(T/V)/Dw")
#        self.phasefigure.set_xlabel("w")
                

class mdOutCpp(DCVizPlotter):
    
    nametag = "ignisPos\d+?\.arma"
    
    transpose = True
    armaBin = True
    
    isFamilyMember = True
    
    loadSequential = True
#    loadLatest = True
    ziggyMagicNumber = 50
    
#    stack = "H"
    
    getNumberForSort = lambda _s, a: int(find("ignisPos(\d+?)\.arma", a)[0])
    
    
    c = [numpy.random.uniform(size=(3,)) for i in range(10)]
    

#    figMap = {"fig": ["subfigure", "eventFigure"]}
    
    
    def plot(self, data):

        x, y = data.data
        

        _colors = ["r", '0.75']
        _sizes = [100, 200]        

        
        _c = [_colors[i%2] for i in range(2)]
        _s = [_sizes[i%2] for i in range(2)]   
        
#        with open("/home/jorgehog/tmp/coolMiddle.arma", 'rb') as f:
#            _data = self.unpackArmaMatBin(f)
#            sliceCorner = [_data[0, 0], _data[0, 1]]
#            sliceSize1 = _data[1, 0]-_data[0, 0] 
#            sliceSize2 = _data[1, 1]-_data[0, 1]
#            self.subfigure.add_patch(pylab.Rectangle(sliceCorner, sliceSize1, sliceSize2, fill=False, color=numpy.random.rand(3,)))
        
        self.subfigure.scatter(x, y, color=_c, s=_s, edgecolor='black')
#        self.subfigure.axes.fill_between([0, Lx], [Ly, Ly], [(1-W)*Ly, (1-W)*Ly], facecolor='red', alpha=0.1)
#        self.subfigure.axes.fill_between([0, Lx], [W*Ly, W*Ly], facecolor='red', alpha=0.1)
#        self.subfigure.axes.fill_between([0, Lx], [(1-W)*Ly/2, (1-W)*Ly/2], [(1+W)*Ly/2, (1+W)*Ly/2], facecolor='blue', alpha=0.1)
        
        __N = self.getNumberForSort(self.filepath)
        self.subfigure.set_title(str(__N))
        self.subfigure.axes.set_xlim([0, 1])
        self.subfigure.axes.set_ylim([0, 1])
        
        
#        try:
#            _data = self.unpackArmaMatBin(open("/home/jorgehog/tmp/mdEventsOut.arma", "rb")).transpose()
#        except:
#            print "failed to load"
#            return
#        _N = re.findall("solver.*?N\s*\=\s*(.+?);", self._cfg, re.DOTALL)[0]
#        T0 = re.findall("mainThermostat.*?bathTemperature\s*\=\s*(.+?);", self._cfg, re.DOTALL)[0]
#        m =  re.findall("Events.*?Thermostats.*?temperatureScaleFactorWarm\s*\=\s*(.+?);", self._cfg, re.DOTALL)[0]
#        
#        T0, _N, m = float(T0), int(_N), float(m)
#        for l, col in enumerate(dataGenerator(_data)):
#            _x = numpy.where(col != 0)
#            self.eventFigure.plot(numpy.arange(__N), col[:__N]/T0, c=self.c[l])
#        self.eventFigure.axes.set_xlim([0, _N-1])
#        self.eventFigure.axes.set_ylim([0, T0*m*1.2])

class KMC_LAMMPS(DCVizPlotter):

    #The name tag of the files to load (of course dumped by lammpswriter.h)    
    nametag = "kMC(\d+)\.lmp"

    #Specify that they are binary
    fileBin = True    
    
    #Specify the sizes of the different header components (LAMMPS)
    binaryHeaderBitSizes = [4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4]
    
    #The number of columns in the binary data file is given as the 11'th header element
    nColsFromHeaderLoc = 11


    #Tells DCViz to load all files in the folder matching the nametag
    isFamilyMember = True
    
    #Instead of loading all family files at once, they will be loaded one by one
    #Can also set this to loadLatest, at which the time stamp of the files will be checked
    #and the last one will be loaded (useful for when running paralell to the simulation)
    loadSequential = True  
    
    #Only reads every 10'th frame.. requested by Sigve, hence the name...
    ziggyMagicNumber = 10 
    
    #Makes font sizes bigger by default (can be further customized)
    hugifyFonts = True

    
    def plot(self, data):
        
        #Potential energy is stored in the 5'th column
        potentialEnergy = data[4]
        print potentialEnergy
        
        #histogram shizz.. self.subfigure and self. figure
        #is the default figure, if more is needed, use this format in the class header:
        #figMap = {'fig1': ['sfig1', 'sfig2', ...], 'fig2': [...], ...}
        self.subfigure.hist(potentialEnergy, bins=40, normed=True)
        
        self.subfigure.set_xbound(0)
        self.subfigure.set_ylim(0, 0.4)
                
        self.subfigure.set_xlabel("Binding energy [E0]")
        self.subfigure.set_ylabel("Density of States")
        

class nucleationHistograms(DCVizPlotter):
    
    nametag = "nucleation(\d+)\.txt"

    transpose = True
    
    hugifyFonts = True
    
    transparent = True
    
    labelSize = 40
    
    def plot(self, data):
        
        print data.shape
        bins, counts = data        
        
        db = bins[1] - bins[0]

        counts /= counts.sum()*db
        
        self.subfigure.bar(bins - db/2, counts, width=db, linewidth=0, facecolor='r')
        
        self.subfigure.set_xbound(0)
        
        ax = self.subfigure.axes
        
        yticks, xticks = pylab.getp(ax, 'xticklabels'), pylab.getp(ax, 'yticklabels')
        
        pylab.setp(yticks, weight="extra bold")
        pylab.setp(xticks, weight="extra bold")
        
#        self.subfigure.set_xlabel("Binding Energy [E0]")
#        self.subfigure.set_ylabel("Density of States")        

class KMC_densities(DCVizPlotter):
    
    nametag = "stateDensity(\d+)\.arma"
    
    armaBin = True

    def plot(self, data):
        
        e, DOS, visit, idx  = data       
        
        print DOS.min(), visit.min(), visit.max(), visit.min()/visit.max()
        
#        DOS *= exp(e)
        
        if DOS.max() != 0:        
            DOS /= DOS.max()
        if visit.max() != 0:
            visit /= 2*visit.max()
        
        self.subfigure.plot(idx, DOS, label="DOS")
        self.subfigure.plot(idx, visit, label="visit")
        self.subfigure.plot([0, idx[-1]], [visit.mean(), visit.mean()], '--r')
        self.subfigure.set_title("flatness = %g" % (visit.min()/visit.mean()))
#        self.subfigure.set_xlim(e.min(), e.max())        
        
        legend()
                

class IGNIS_EVENTS(DCVizPlotter):
    
    nametag = "ignisEventsOut*"
    
    armaBin = True

    def plot(self, data):
        
        T, N, E, t = data
        
        i = where(t!=0)
        
        print T[0]
        self.subfigure.plot(T[i])
        legend()

"""
    figMap = {"nfig" : ["N", "T", "E"]} 
    
    stack = "V"    
    
    def plot(self, data):
        
        T, N, E = data

        T = T[numpy.where(T != 0)]
        N = N[numpy.where(N != 0)]     
        E = E[numpy.where(E != 0)]

        s = min([len(N), len(T), len(E)])        
        
        N = N[:s]
        T = T[:s]        
        E = E[:s]
        
        cN = numpy.cumsum(N)/numpy.cumsum(numpy.ones(len(N)))
#
#        r = 100000
#        cN = numpy.zeros(len(N)-2*r)
#        for i in range(r, len(N) - r):
#            cN[i-r] = numpy.mean(N[i-r:i+r])


        dN = cN[1:] - cN[:-1]
        dT = T[1:] - T[:-1]        
        
        cut = 10000
        dN = dN[cut:]
        dT = dT[cut:]        
        
#        dNdT = dN/dT        
        
        self.N.plot(N)
        self.N.set_ylabel(r"avg N")
        self.N.set_xlabel("t")
        self.N.set_xbound(0)
        self.N.set_ybound(0)
        self.N.set_title("T = %g" % T[-1])
"""
  
import glob
class molsim2(DCVizPlotter):
    
    nametag = "ChemPot|AvPress"
    
    isFamilyMember = True
    
    def plot(self, data):
        
        rho = "0.001 0.003 0.006 0.009 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.33 0.5 0.6 0.7 0.8 0.9"
        rho = [float(r) for r in rho.split()]
        
        s = 0
        e = len(rho)-1        
        
        chem, pres = data
        
        self.subfigure.plot(chem.data[s:e], pres.data[s:e])
        
        for c, p, r in zip(chem.data[s:e], pres.data[s:e], rho[s:e]):
            
            self.subfigure.text(c, p, str(r))

class molsim3(DCVizPlotter):
    
    nametag="lj\.densplot"

#    isFamilyMember = True

    def plot(self, data):
        
#        self.dictify(data)
        
        n, d1, d2 = data
        
        self.subfigure.plot(n, d1)
        self.subfigure.plot(n, d2)
    
    
class MD_OUT(DCVizPlotter):
    
    nametag = "MD_out\d*\.dat"

    skipRows = 2
    isFamilyMember = True
    
#    loadLatest = True    
    loadSequential = True
        
    
    getNumberForSort = lambda _s, a : int(find("MD\_out(\d+)\.dat", a)[0])  
    
    colorTable = ['b', 'r', 'g', 'k', '0.5']
    
    figMap = {"fig": ["subfigure", "veldist"]}
    
    stack = "H"
    
    def plot(self, data):
      
        lx, ly, T = [float(x) for x in self.skippedRows[0].split()]
        Tz = self.skippedRows[1]
        
        
        
        C, X, Y, Vx, Vy = data
        
        for c, x, y in zip(C, X, Y):
            self.subfigure.plot(x, y, self.colorTable[int(c)] + "o")
        
        V = numpy.sqrt(Vx**2 + Vy**2)
        
        try:
            self.veldist.hist(V[numpy.where(V > 0)], bins=20, facecolor="green")
        except:
            pass
        
        
        self.subfigure.axes.set_xlim([0, lx])
        self.subfigure.axes.set_ylim([0, ly])
        self.subfigure.set_xlabel("x")
        self.subfigure.set_ylabel("y")
        self.veldist.set_xlabel("hist(v)")
        self.subfigure.set_title(os.path.split(self.filepath)[-1].split(".")[0].split("_")[1].replace("out", ""))
        self.veldist.set_title("T = %g / %.2f %.2f %.2f" %  ((T,) + tuple(eval(Tz))))
        
    

class EnergyTrail(DCVizPlotter):
    
    nametag = "\w+.trailingE.arma"
    
    armaBin=True
    isFamilyMember=True
    
    figMap = {"fig1": ['eFig']}
    c = ["#C0C0C0", '#008000', '#008000'] 
    l = ["-", "--", "-"]
    
    def plot(self, data):

        eDMC = 65.7
        
        for i, data_ in enumerate(data):
            name = self.familyFileNames[i].split(".")[0].replace("_", " ")[1:]
            self.eFig.plot(numpy.cumsum(data_.data)/numpy.linspace(1,data_.m,data_.m),
                           self.l[i], c=self.c[i], label=name)
                           
        self.eFig.plot([1,data_.m], [eDMC, eDMC], 'k-.', label="DMC Energy")
        self.eFig.legend()
        
        self.eFig.set_xlabel('Cycle')
        self.eFig.set_ylabel(r'$\langle E\rangle $', rotation=0)
        
        self.eFig.axes.get_yaxis().get_label().set_fontsize(20)
        self.eFig.axes.get_xaxis().get_label().set_fontsize(20)

        formatter = ticker.ScalarFormatter(useOffset=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-3, 2))
        self.eFig.axes.get_xaxis().set_major_formatter(formatter)


class Blocking(DCVizPlotter):
    
    nametag = "blocking_\w+_out\d*\.dat"
    figMap = {"Fig": ["blockFig"]}
    
    nameMap = {"0": r"$\alpha$", "1": r"$\beta$", "": ""}
    
    def plot(self, data):
        
        Fig, blockFig = self.Fig, self.blockFig
        blockSize, error = data
        
#        fileName = os.path.basename(self.filepath)
#        title = "Blocking data %s" % (fileName.split("_")[1] + " %s" % \
#              (self.nameMap[re.findall("blocking_\w+_out(\d*)\.dat", self.filepath)[0]]))
        blockFig.plot(blockSize, error, 
                      '*', color='#008000', markersize=10)  
        
        formatter = ticker.ScalarFormatter(useOffset=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-3, 4))
        blockFig.axes.get_yaxis().set_major_formatter(formatter)
        blockFig.axes.get_yaxis().get_label().set_fontsize(30)
        blockFig.axes.get_yaxis().offsetText.set_size(20)
        blockFig.axes.get_xaxis().get_label().set_fontsize(20)

        
#        blockFig.set_title(title)
        blockFig.set_xlabel(r'Block size')
        blockFig.set_ylabel(r'$\sigma$', rotation=0)
        

class DMC_OUT(DCVizPlotter):
    
    nametag = "DMC_out\.dat"
    figMap = {"Fig": ["N_plot"], "Fig2": ["E_plot"]}
    dt = 0.001
        
        
    def plot(self, data):
        
        E, Eavg, N, Navg, ET = data
        
        t = numpy.linspace(0, self.dt*(len(E) - 1), len(E))
        N_plot, E_plot = self.N_plot, self.E_plot
        
        lw=2
        
        # E PLOTS
        E_plot.plot(t, E, 'k--', label="dmc E", linewidth=lw)
        E_plot.plot(t, ET, color="#008000", label="trial E", aa=True
        )
        
        E_plot.legend()
#        E_plot.set_title('Energy convergeance')
        E_plot.set_xlabel(r'$\tau = n*\delta \tau$ [s]')
        E_plot.set_ylabel(r'E [Ha]')
        E_plot.ticklabel_format(useOffset=False, axis='y')
        
        E_plot.axes.get_yaxis().get_label().set_fontsize(20)
        E_plot.axes.get_xaxis().get_label().set_fontsize(20)
        
        # N PLOTS
        N_plot.plot(t, N, '#008000')
        N_plot.plot(t, Navg, 'k--', linewidth=lw)
        
        N_plot.set_ylabel(r"N_W($\tau$)")
        N_plot.axes.get_xaxis().set_visible(False)
#        N_plot.set_title('Walker population')
        N_plot.ticklabel_format(useOffset=False, axis='y')
        
        N_plot.axes.get_yaxis().get_label().set_fontsize(20)

class radial_out(DCVizPlotter):
    
    nametag = "radial_out.+\.arma"
    figMap = {"fig1":["radialFig"]}
    
    armaBin = True
    isFamilyMember=True
    familyName = "radial dist"
    
    def plot(self, data):
        cut=4
        xScale = 1
        
        piFac = 2*3.141592;
        for d in data:
            d.data *= piFac
   
        silent = True
        
        color = ['#008000', "0.5", "k", '#008000']
        style = ['-', '-.', '--', '.']
        max_edge = 0   
        maxCut = 0
        j = 0
        
        path, name = os.path.split(self.filepath)
        
        try:
            yFile = open(os.path.join(path, "yMIN.dat"), 'r')
            yMax = float(yFile.read())
            print yMax
        except:
            yMax = None
        
        vmc=None
        dmc=None
        
        pureOnly = False
        superPose = False  
        dmcEdge = None
        vmcEdge = None
        
        r2 = False
        if "Atoms" in self.familyFileNames[0]:
            r2 = True
        
        if superPose:
            yMax = 1.2

        for i in range(len(data)):
            
            if superPose:
                edge = 1
            else:
                edge = float(re.findall("edge(\d+\.?\d*)\.arma", self.familyFileNames[i])[0])
                if not silent: print "n_p= ", data[i].data.sum()*edge/(data[i].n-1), "?"
            
            if edge > max_edge:
                max_edge = edge;
            
            if "dmc" in self.familyFileNames[i]:
                method = "dmc"
                dmc = data[i].data
                
                dmcEdge = edge;

                if superPose:
                    dmc = dmc/dmc[cut:].max()
                
                last = dmc
                
            elif "vmc" in self.familyFileNames[i]:
                method = "vmc"
                vmc = data[i].data
                vmcEdge = edge
                
          
                if superPose:
                    vmc = vmc/vmc[cut:].max()
     
                last = vmc
                
                
            r = numpy.linspace(0, edge, data[i].n)
            
            if r[cut] > maxCut:
                maxCut = r[cut]            
            
            if not pureOnly:
                if "QDots3D" in self.familyFileNames[i] or "Diatom" in self.familyFileNames[i]:
                    r.resize((data[i].n, 1))
                    last[1:] /= r[1:]**2
                    last[0] = last[1]
                    
                    if "vmc" in self.familyFileNames[i]:
                        vmc = last
                    else:
                        dmc = last
           
                self.radialFig.plot(r, last, style[i%2], label=method.upper(), color=color[i%2]);
           
            if vmc is not None and dmc is not None:
                if vmcEdge != dmcEdge:
                    print "Warning. Ploting pure dist for mismatching edges. %f != %f" % (vmcEdge, dmcEdge)
                else:
                    print "Pure success"
                    
                if pureOnly:
                    pureC = color[j%len(color)]
                    pureS = style[j%len(color)]
                    j += 1
                    pLabel=None

                    if superPose:
                        pLabel = re.findall("out_(.+?)[vd]mc", self.familyFileNames[i])[0]
                        pLabel = re.sub("(\d)c(\d)", "\g<1> \g<2>", pLabel)
                        pLabel = re.sub("QDots\d+", "", pLabel)
                    
                else:
                    pureC = 'k'
                    pureS = "--"
                    pLabel= "Pure"
                    
                pure = 2*dmc - vmc

                if superPose:
                    pure = pure/pure[cut:].max()
                
                self.radialFig.plot(r, pure, pureS, color=pureC, label=pLabel)
                vmc = None
                dmc = None
           
        
#        if not pureOnly:
        self.radialFig.legend()    
        self.radialFig.axes.set_xlim(maxCut, xScale*max_edge)
#        self.radialFig.axes.set_xlim(maxCut, 5)
        self.radialFig.set_xlabel('r')
        
        if r2:
            self.radialFig.set_ylabel(r'$r^2\rho(r)$')
        else:
            self.radialFig.set_ylabel(r'$\rho(r)$')

#        locator = self.radialFig.axes.get_yaxis().get_major_locator()
#        self.radialFig.axes.set_ylim(locator.autoscale()/2)
        if yMax is not None:
            self.radialFig.set_ylim(0, yMax)
        self.radialFig.axes.set_ybound(0)
        
        self.radialFig.axes.get_yaxis().get_label().set_fontsize(30)
        self.radialFig.axes.get_xaxis().get_label().set_fontsize(30)
        
        

class dist_out(DCVizPlotter):
    
    nametag = "dist_out.+\.arma$"
    
    figMap = {"fig1": [], "fig2": ["subfig1"]}

    armaBin = True
    isFamilyMember = True
    familyName = "dist 2D"        
        
    stack = "H"
    
    dmcOnly = True
    vmcOnly = False
    
        
    def plot(self, data):
    
        silent = False
          
        edge = float(re.findall("_edge(.+?)\.arma", self.familyFileNames[0])[0])    
    
        if len(data) == 1:
            if not silent: print "length 1 data"
            dist = data[0].data
        elif len(data) == 2:
            if not silent: print "len2 data"
            edge_2 = float(re.findall("_edge(.+?)\.arma", self.familyFileNames[1])[0])
        
            if not silent:
                if edge != edge_2:
                    print "Bin edges does not match. %s != %s" % (edge, edge_2)

            for i in range(2):
                if "vmc" in self.familyFileNames[i]:
                    vmcDist = data[i].data
                elif "dmc" in self.familyFileNames[i]:
                    dmcDist = data[i].data
            
        
            if self.dmcOnly:
                try:
                    dist = dmcDist
                except: 
                    if not silent: print "\n\nWarning: No DMC data found. Attempting to load VMC data.\n\n"
                   
                    dist = vmcDist
                   
            elif self.vmcOnly:
                dist = vmcDist
            else:
                try:
                    dist = 2*dmcDist - vmcDist
                    print "pure success!", dmcDist.sum()*(edge/100)**2, vmcDist.sum()*(edge/100)**2 
                except:
                    raise Exception("Supplied dist files does not match a VMC+DMC pair:  \n%s \n%s" % (self.familyFileNames[0], self.familyFileNames[1]))
        else:
            raise Exception("More than two distributions loaded in given folder")
        
        
        origLen = len(dist)
        distMid = dist[:, origLen/2]
        crit = numpy.where(distMid > 0.1*distMid.max())[0]

        x, y = numpy.meshgrid(crit, crit)
        dist = dist[x, y]
        
        ax = Axes3D(self.fig1)#, self.subfig0.get_position())
        
        r = numpy.linspace(-edge, edge, origLen)
        
        X, Y = numpy.meshgrid(r, r)
        X = X[x, y]
        Y = Y[x, y]
        C = cm.Greens
        
        if dist[len(dist)/2:, :].sum() < dist[0:len(dist)/2, :].sum():
            print "flipping"
            print dist[:, 0].sum()
            dist = numpy.flipud(dist)
            print dist[:, 0].sum()
        ax.plot_surface(X, Y, dist, rstride=1, cstride=1, cmap=C, linewidth=0)


        if "DoubleWell" in self.familyFileNames[0]:
            zdir = "y"
            offset=-ax.get_ylim()[0]
        else:
            zdir = "x"
            offset = ax.get_xlim()[0]
        
        cset = ax.contour(X, Y, dist, zdir=zdir, offset=offset*1.05, color='#008000', levels=[0]) 
           
        ax.set_zlim(0, dist.max())

        ax.set_ylabel("x")
        ax.set_xlabel("y")
        ax.view_init(30, -65)
#        ax.view_init(0, 90)
        
#        print dist.shape
#        dist = zoom(dist, 3)
#        extent = [-newEdge, newEdge, -newEdge, newEdge]
        self.subfig1.axes.contourf(X, Y, dist, zdir='z', cmap=C)
        self.subfig1.set_xlabel("x")
        self.subfig1.set_ylabel("y", rotation=0)
        self.subfig1.axes.get_yaxis().get_label().set_fontsize(20)
        self.subfig1.axes.get_xaxis().get_label().set_fontsize(20)
        
#        self.subfig0.axes.get_xaxis().set_visible(False)
#        self.subfig0.axes.get_yaxis().set_visible(False)
        
        
#        self.subfigHist3D.set_ylabel(r'y')
#        self.subfigHist3D.set_xlabel(r'x')
#        self.fig1.colorbar(im)

#        self.subfigDist1d.set_ylabel(r'$|P(r)|^2$')


class R_vs_E(DCVizPlotter):
    
    nametag = "R\_vs\_E.*?\.dat"
    
#    figMap = {"fig":["sfigV", "sfigE"]}
    figMap = {"fig":["sfigE"]}
    stack="H"

    def plot(self, data):
        R, Ep, Ec, Ek = data
        
        R, Ep, Ec, Ek = zip(*sorted(zip(R, Ep, Ec, Ek), key=lambda x: x[0]))
      
        
#        E = Ep + Ec + Ek
        
        R = numpy.array(R)
#        print R.shape        
#        print E[-1]
#        i = numpy.where(abs(E) > 2*abs(E[-1]))[0][-1]
        i = numpy.where(R > 0.5)[0][0]        
        print i
#        R = R[i:]
        Ep = numpy.array(Ep)
        Ec = numpy.array(Ec)
        Ek = numpy.array(Ek)
        
#        self.sfigE.plot(R[i:], (Ep + Ek + Ec)[i:], '*', color='#008000', label=r"$\langle E\rangle$")
#        self.sfigE.set_xlabel("R")
#        self.sfigE.legend()
        
        LJ = numpy.vectorize(lambda r, eps, sig: 4*eps*(sig**12/r**12 - sig**6/r**6))
        
        
        self.sfigE.plot(R, LJ(R, 1., 1.), color='#008000', label=r"$V(R)$")
        self.sfigE.set_xlabel("R")
        self.sfigE.set_ylim(-1.5, 2)
        self.sfigE.axes.set_xbound(0.5)
        self.sfigE.legend()        
        
#        self.sfigV.plot(R, Ep + Ec, '*', color='#008000', label=r"$\langle V\rangle$")
#        self.sfigV.set_xlabel("R")
#        self.sfigV.set_ylim([-3, -1])
#        self.sfigV.legend()
#                
        
#        self.sfigV.axes.set_xbound(0.25)
#        self.sfigE.axes.set_xbound(0.25)
#        self.sfigV.axes.get_xaxis().set_visible(False)

class Scaling(DCVizPlotter):
    figMap = {"f1":["fig", "fig2"], "f2":["fig3", "fig4"]}
    nametag = "scalePlot\.dat"
    
    stack = "H"
    
    def plot(self, data):
        
#        Q2D_N = numpy.array([2, 6, 12, 20, 30, 42, 56])
#        Q2D_T = numpy.array([0.2, 1.06, 4.13, 12.71, 33.9, 50.76, 107])
        Q2D_N = numpy.array([2, 6, 12, 20, 30, 42, 56])
        Q2D_T = numpy.array([0.17, 1.05, 4.42, 13.25, 34.63, 80.24, 165.77])
        
#        Q3D_N = numpy.array([2, 8, 20])
#        Q3D_T = numpy.array([0.38, 2.42, 15.3])
        Q3D_N = numpy.array([2, 8, 20])
        Q3D_T = numpy.array([0.38, 2.42, 16.59])
        
#        Atoms_N = numpy.array([2, 4, 10, 12, 18, 36])
#        Atoms_T = numpy.array([0.4, 0.88, 3.92, 5, 13.3, 42.3])
        Atoms_N = numpy.array([2, 4, 10, 12, 18, 36])
        Atoms_T = numpy.array([0.42, 0.85, 3.94, 5.67, 13.5, 66])
        
#        Molecules_N = numpy.array([2, 6, 8, 10, 12, 14, 16])
#        Molecules_T = numpy.array([0.47, 2.02, 3.34, 5, 7.33, 10.43, 13.4])
        Molecules_N = numpy.array([2, 6, 8, 10, 12, 14, 16])
        Molecules_T = numpy.array([0.49, 2.02, 3.27, 5, 7, 9.42, 12.24])
        
        
        from scipy.stats import linregress as l2
          
            
              
        
        toPlot = [[Q2D_N, Q2D_T], [Q3D_N, Q3D_T], [Atoms_N, Atoms_T], [Molecules_N, Molecules_T]]
        l = ["Qdots 2D", "Qdots 3D", "Atoms", "Molecules"]
        c = ['#008000', 'k', '0.5', "#008000"]
        style = ["-*", "-^", "-o", "-d"]        
        
        for i, NT in enumerate(toPlot):
            N, T = NT
            k = numpy.where(N <= 20)
            
            if len(N[k]) != 3:
                slope, intercept, r_value, p_value, std_err = l2(numpy.log(N[k][1:]), numpy.log(T[k][1:])) 
            else:
                r_value = 1
                slope = (numpy.log(T[k][2]) - numpy.log(T[k][1]))/(numpy.log(N[k][2]) - numpy.log(N[k][1]))
            print slope, r_value, l[i]
            
            self.fig.loglog(N[k], T[k], style[i], color=c[i], label=l[i])     
            self.fig3.plot(N[k], T[k], style[i], color=c[i], label=l[i])
        self.fig.legend(loc=2)
        self.fig.set_ylabel("t[s]")
        self.fig.set_xlabel("N")
        self.fig3.legend(loc=2)
        self.fig3.set_ylabel("t[s]")
        self.fig3.set_xlabel("N")
        
        self.fig2.loglog(Q2D_N, Q2D_T, "-*", color=c[0], label=l[0])
        self.fig2.loglog(Atoms_N, Atoms_T, "-o", color=c[1], label=l[2])
        self.fig2.legend(loc=2)
        self.fig2.set_ylabel("t[s]")
        self.fig2.set_xlabel("N")
        
        self.fig4.plot(Q2D_N, Q2D_T, "-*", color=c[0], label=l[0])
        self.fig4.plot(Atoms_N, Atoms_T, "-o", color=c[1], label=l[2])       
        self.fig4.legend(loc=2)
        self.fig4.set_ylabel("t[s]")
        self.fig4.set_xlabel("N")
        
        
class E_vs_w(DCVizPlotter):
    
    stack = "H"    
    
    nametag = "E\_vs\_w\.dat"
    figMap = {
              "f2":["s2", "se2"],
              "g2" : ["sg2"],
              "f6":["s6", "se6"],
              "g6" : ["sg6"],
              "f12":["s12", "se12"],
              "g12" : ["sg12"], 
              "f20":["s20", "se20"], 
              "g20" : ["sg20"],            
              "f30":["s30", "se30"], 
              "g30" : ["sg30"],            
              "f42":["s42", "se42"],
              "g42": ["sg42"]
              }
              
#    figMap = {"f1": ["sg"]}
 
    def plot(self, data):
        
        N, W, E, E_K, E_O, E_C = data
        
        I = {"2": [-1, -1], "6": [-1, -1], "12": [-1, -1], "20": [-1, -1], "30": [-1, -1], "42": [-1, -1]}
        
        k = 0
        for n_p in N:
            i = str(int(n_p))
            if I[i][0] == -1:
                I[i][0] = k
           
            I[i][1] = k
            
            k+=1
         
         
        for N in ["2", "6", "12", "20", "30", "42"]:
            i1 = I[N][0]
            i2 = I[N][1]+1

            w = W[i1:i2][numpy.where(E[i1:i2] < 1000)]
            e = E[i1:i2][numpy.where(E[i1:i2] < 1000)]
            ek =E_K[i1:i2][numpy.where(E[i1:i2] < 1000)]
            eo = E_O[i1:i2][numpy.where(E[i1:i2] < 1000)]
            ec = E_C[i1:i2][numpy.where(E[i1:i2] < 1000)]

            subfig = eval("self.s%s" % N)
            subfig2 = eval("self.se%s" % N)
            subfigG = eval("self.sg%s" % N)
#            subfigG = eval("self.sg%s" % "")
#            
 

            subfig.plot(w, ec/e, ".", color='#008000', label="Ecol/E")
            subfig.plot(w, eo/e, "^", color='#008000', label="Eosc/E")
            subfig.plot(w, ek/e, "*", color='#008000', label="Ekin/E")
#            subfig.set_title("N = %s" % N)     
            subfig.legend()
            subfig.set_ylim([0, 1])
            subfig.set_xlabel("$\omega$")
            
            subfig2.plot(w, e/w, "+", color='#008000', label="E/$\omega$")
            subfig2.plot(w, ec/w, ".", color='#008000', label="Ecol/$\omega$")
            subfig2.plot(w, eo/w, "^", color='#008000', label="Eosc/$\omega$")
            subfig2.plot(w, ek/w, "*", color='#008000', label="Ekin/$\omega$")
         
            subfig2.legend()
            subfig2.axes.set_ybound(0)
            subfig2.set_xlabel("$\omega$")

            Vtot = (eo + ec)
         
#            Vtot /= Vtot.mean()

#            ek /= ek.mean()
            np = int(N)
            
            alpha = 1.3
            beta = 1.83
            
            ek /= np**alpha
            Vtot /= np**beta
            
            subfigG.plot(Vtot, ek, "o", color="0.5", mfc="None", label="QMC Result")
            subfigG.set_ylabel("Ekin/N^%.2f" % alpha)
            subfigG.set_xlabel("(Eosc + Ecol)/N^%.2f" % beta)
#            subfigG.set_title("N = %s" % N)
                        
            from scipy.stats import linregress as l2
            n = len(ek)
            
            treshHigh = n/2
            treshLow = n/4
            
            slope, intercept, r_value, p_value, std_err = l2(Vtot[treshHigh:], ek[treshHigh:])
            lFit = intercept + slope*Vtot

            subfigG.plot(Vtot, lFit, '-', color='#008000', label=r"r2 = %g. a = %f" % (r_value**2, slope))
            
            slope, intercept, r_value, p_value, std_err = l2(Vtot[:treshLow], ek[:treshLow]) 
            lFit = intercept + slope*Vtot

            subfigG.plot(Vtot, lFit, 'k--', label=r"r2 = %g. a = %f" % (r_value**2, slope), linewidth=2)
            
            subfigG.legend(loc=2)
            subfigG.axes.set_ybound(0)
            subfigG.axes.set_xbound(0)
    
        
class testBinFile(DCVizPlotter):
    
    nametag = "testBin.+\.arma"
    figMap = {"fig": ["subfig"]}

    fileBin = True
    Ncols = 2
    
    skipRows = 2    

    def plot(self, data):
        x, y = data
        self.subfig.plot(x, y, '*')
        
class MIN_OUT(DCVizPlotter):
    
    nametag = "ASGD_out\.dat"
    figMap = {"E_fig"    : ["E_plot"], 
              "step_fig" : ["step_plot"],
              "param_fig": ["param_plot", "grad_plot"]}
    
    indexmap = {0: r"\alpha", 1: r"\beta"}
    l = ["-", "--"]
    c = ['#008000', 'k']  
    
    transpose = True
    
    def plot(self, data):

        n_params = (self.Ncols - 3)/2


        E, Eavg, step = dataGenerator(data[:3])
        
        E_plot, step_plot, param_plot, grad_plot = self.E_plot, self.step_plot, \
                                                    self.param_plot, self.grad_plot
        
        #~ E PLOTS
        E_plot.plot(E, self.l[0], c=self.c[0], label="E")
        E_plot.plot(Eavg, self.l[1], c=self.c[1], label="average E")
        
        E_plot.legend()
        #E_plot.set_title('Energy convergeance')
        E_plot.set_xlabel(r'Cycle')
        E_plot.set_ylabel(r'E [Ha]')
        E_plot.axes.get_yaxis().get_label().set_fontsize(30)
        E_plot.axes.get_xaxis().get_label().set_fontsize(30)
        E_plot.ticklabel_format(useOffset=False, axis='y')
        
        #~ Step plot
        step_plot.plot(step, self.l[0], c=self.c[0])
        #step_plot.set_title('Step length')
        step_plot.set_xlabel('Cycle')
        step_plot.set_ylabel('Step')
        step_plot.axes.get_yaxis().get_label().set_fontsize(30)
        step_plot.axes.get_xaxis().get_label().set_fontsize(30)

        lw = 2
        #~ Param plots
        for i in range(0, 2*n_params, 2):
            param_plot.plot(data[3 + i], self.l[i/2], c=self.c[i/2],label=r'$%s$' % self.indexmap[i/2], linewidth=lw)    
            grad_plot.plot(data[4 + i], self.l[i/2], c=self.c[i/2], linewidth=lw)
        
        param_plot.set_ylabel(r'$\alpha_i$')
        #param_plot.set_title('Variational parameters')
        param_plot.axes.get_xaxis().set_visible(False)
        
        if n_params > 1:
            param_plot.legend()
        param_plot.axes.get_yaxis().get_label().set_fontsize(30)
        param_plot.axes.get_xaxis().get_label().set_fontsize(30)
        param_plot.axes.get_yaxis().labelpad = 20
        
        #grad_plot.set_title('Energy derivatives')
        grad_plot.set_ylabel(r'$\frac{\partial E}{\partial \alpha_i}$', rotation=0)
        grad_plot.set_xlabel('Cycle')
        grad_plot.axes.get_yaxis().get_label().set_fontsize(30)
        grad_plot.axes.get_xaxis().get_label().set_fontsize(30)
            

    
        
if __name__ == "__main__":
    print "Invalid usage: Use the DCVizWrapper for terminal usage."
