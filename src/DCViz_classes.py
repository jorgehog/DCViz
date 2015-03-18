
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

try:
    from scipy.stats import linregress
except:
    pass

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

class hist_stuff(DCVizPlotter):

    nametag = "hist\.arma"

    armaBin = True

    def plot(self, data):
        self.subfigure.plot(data.data)

class acf_height(DCVizPlotter):

    nametag = "acf\.arma"

    armaBin = True

    def plot(self, data):

        self.subfigure.plot(data.data/data.data.max())
        self.subfigure.set_xlim(0, 30)
        self.subfigure.set_xlabel("dx")
        self.subfigure.set_ylabel("acf(dx)")

class SolidOnSolidHeights(DCVizPlotter):

    nametag = "SOSHeights\.arma"

    armaBin = True

    figMap = {"fig3D" : []}

    def plot(self, data):
        ax = Axes3D(self.fig3D)

        hist = data.data - data.data.min() + 1

        elements = data.n*data.m

        xpos, ypos = np.meshgrid(arange(data.n), arange(data.m))

        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(elements)
        dx = np.ones_like(zpos)
        dy = dx.copy()
        dz = hist.flatten()

        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, linewidth=0)



class heighMap(DCVizPlotter):

    nametag = "heights\.arma"

    armaBin = True

    def plot(self, data):

        fs = rcParams["figure.figsize"]
        ar = fs[1]/fs[0]
        y_max = int(ar*data.n)

        Z = data.data - data.data.min()

        if (y_max/2 > Z.max()):
            Z += y_max/2;

            for i, zi in enumerate(Z):
                self.subfigure.scatter(i+0.5, zi, marker='s', color='r', s = 0.09*(data.n/fs[0])**2)

                self.subfigure.set_ylim(0, y_max)
        else:
            self.subfigure.set_ybound(0)

        self.subfigure.bar(arange(data.n), Z, width=1.0, linewidth=0)

        self.subfigure.set_xlim(0, data.n-1)

class testStuff(DCVizPlotter):
    
    nametag = "heights\.arma"
    
    armaBin = True
    
    def plot(self, data):

        self.subfigure.plot(data.data, 'b*')


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
    
    figMap = {"fig" : ("bohrfigure"), "fig2": "rfigure", "fig3": "rrelfigure"}
    
    transpose = True    
    
    def plot(self, data):
        
        for i, _data in enumerate(data):
            N = re.findall(self.nametag, self.familyFileNames[i])[0]
            
            w, alpha, beta, E_v, T_v, vho_v, vcol_v, r_v, rij_v, E_d, T_d, vho_d, vcol_d, r_d, rij_d = _data
        
            V_v = vho_v + vcol_v
            V_d = vho_d + vcol_d
            
            l = "N=%s" % N              
            
            print l
            for o, E in zip(w, T_d+V_d):
                print "%g   %g"  % (o, E)
            print "------------------------------"
            
            GAMMA_v = T_v/V_v
            GAMMA_d = T_d/V_d
                       
            self.bohrfigure.plot(w, GAMMA_d, label=l + "dmc", marker='^', linestyle="--")
            self.bohrfigure.plot(w, GAMMA_v, label=l + "vmc", marker='^', linestyle="--")
            self.rfigure.plot(w, r_d/r_d.max(), "--^", color="k", label="r/%g" % r_d.max())
            self.rfigure.plot(w, rij_d/rij_d.max(), "-*", color="r", label="rij/%g" % rij_d.max())

        self.bohrfigure.legend(loc=4)
        self.bohrfigure.set_xlabel("$\omega$")
        self.bohrfigure.set_ylabel("T/V")

        self.rfigure.legend(loc=1)
        self.rfigure.set_xlabel("$\omega$")
        self.rfigure.set_ylabel("r")
        
        self.rrelfigure.set_xlabel("$\omega$")                
        self.rrelfigure.set_ylabel("r_rel")

        
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


class IGNIS_TEST(DCVizPlotter):

    nametag = "ignis_test\.ign"

    fileBin = True

    binaryHeaderBitSizes = [4, 4]
    nColsFromHeaderLoc = 1

    def plot(self, data):

        for i in range(self.Ncols):
            self.subfigure.plot(data[i], label="%d" % i)

        legend()


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

class quick_conc(DCVizPlotter):

    nametag = "conc(\d+).arma"

    #isFamilyMember = True

    #loadLatest = True

    armaBin = True

    def plot(self, data):

        self.subfigure.plot(data.data)
        self.subfigure.set_ybound(0)

class quick_hist(DCVizPlotter):
    
    nametag = "hist\.npy"    
    
    numpyBin = True
    figMap = {"fig" : ["subfigure", "sf2"]}

    def plot(self, data):
        
        self.subfigure.plot(data.data)
        self.subfigure.set_ybound(0)
        self.sf2.plot(data.data)

class WLMC_C(DCVizPlotter):

    nametag = "wl(\d+)\.dat"
    isFamilyMember = True
    loadSequential = True

    transpose = True

    def plot(self, data):

        print data.shape
        bin, E, Eavg, H, lng = data

        self.subfigure.plot(E, lng)
        self.subfigure.set_title(self.filename)

class KMC_densities(DCVizPlotter):
    
    nametag = "stateDensity(\d+)\.arma"
    
    armaBin = True

    plotvisit = True

    hugifyFonts = True

    def mapIndex(self, indices, index):
        hits = where(indices == index)[0]

        if len(hits) == 0:
    
            minError = 100000000;
            winner = None
            
            for i, idx in enumerate(indices):
                
                if abs(idx - index) < minError:
                    minError = abs(idx - index)
                    winner = i
                
            return winner        
        
        return hits[0];

    def hasAwesome(self):

        if not "canonical_logMTVec.arma" in os.listdir("/tmp/"):
            return False
        return True

    def loadaws(self):

        with open("/tmp/canonical_logMTVec.arma", 'rb') as f:
            return self.unpackArmaMatBin(f)[1:]


    def plot(self, data):
        
        e, DOS, visit, idx  = data       

        # DOS -= DOS.min()
        #
        # if DOS.max() != 0:
        #     DOS /= DOS.max()
        if visit.max() != 0:

            visit *= (DOS.max()/(2*visit.max()))



        if self.hasAwesome():

            aws = self.loadaws()

            self.subfigure.plot(aws, 'r^', fillstyle="none")

            X = linspace(0, len(aws), len(DOS))
        else:
            X = idx

        self.subfigure.plot(X, DOS, 'bs', fillstyle="none")

        self.subfigure.set_xlim(0, X[-1])

   #     self.subfigure.plot(idx[1:], idx[1:]**-0.5*DOS[1]/idx[1])

        if self.plotvisit:
            self.subfigure.plot(X, visit, 'g', label="visit")

            self.subfigure.set_title("flatness = %g" % (visit.min()/visit.mean()))
        self.subfigure.set_xbound(-0.5)
        self.subfigure.set_ybound(-0.5)

        self.subfigure.set_xlabel("s")
        self.subfigure.set_ylabel("$\log(\Omega (s))$")


        try:
            with open(os.path.join(self.path, 'flatness.txt'), 'r') as f:
    
                for line in f:

                    l, u = [int(x) for x in line.split()]
    
                    l = self.mapIndex(idx, l)
                    u = self.mapIndex(idx, u-1) + 1

                    el = X[l]
                    eu = X[u-1] + 1

                    m = visit[l:u].mean()

                    self.subfigure.plot([el, eu], [m, m], 'k--*')
                    
    
            with open(os.path.join(self.path, 'roughness.txt'), 'r') as f:
    
                for line in f:
    
                    l, u = [int(x) for x in line.split()]
    
                    l = self.mapIndex(idx, l)
                    u = self.mapIndex(idx, u - 1) + 1

                    el = X[l]
                    eu = X[u-1] + 1

                    self.subfigure.plot([el, eu], [1./2, 1./2], 'r--*')
        except:
            pass
#        self.subfigure.set_xlim(e.min(), e.max())        
        
        legend(loc=0)
                

class logDOS_test(DCVizPlotter):

    nametag = "canonical_logMTVec\.arma"

    armaBin = True

    def plot(self, data):
        s = arange(len(data))
        self.subfigure.plot(data.data)


class QuasiLoadedIgnis20(DCVizPlotter):

    nametag = "ignisSOS.ign"

    fileBin = True

    binaryHeaderBitSizes = [4, 4]
    nColsFromHeaderLoc = 1

    figMap = {"fig1" : ["sfig", "hfig"]}

    def plot(self, data):
        s, h = data

        self.sfig.plot(s)
        self.hfig.plot(h)

class QuasiLoadedIgnis(DCVizPlotter):

    nametag = "ignisQuasi2Dloaded.ign"

    fileBin = True

    binaryHeaderBitSizes = [4, 4]
    nColsFromHeaderLoc = 1

    figMap = {"fig" : ["subfigure", "subfigure2", "subfigure3"], "fig2" : ["subfigure4"]}


    def plot(self, data):

        S = 100

        if data.m == 8:
            t, h, dh, surfsize, s, hw, eqC, mc = data

            eqC = eqC[::S]
            mc = mc[::S]

        elif data.m == 7:
            t, h, dh, surfsize, s, eqC, mc = data
            hw = h

        elif data.m == 5:
            t, h, dh, surfsize, s = data
            eqC = empty(0)
            mc = empty(0)
            hw = h

        else:
            t, h, dh, surfsize, s, hw = data
            eqC = empty(0)
            mc = empty(0)

        t = t[::S]
        h = h[::S]
        dh = dh[::S]
        hw = hw[::S]


        self.subfigure2.plot(t, h, 'b')
        self.subfigure2.set_xlabel("t")
        self.subfigure2.set_ylabel("h")

        nonZeroDh = where(dh != 0)

        if not nonZeroDh:
            return

        dh = dh[nonZeroDh]
        self.subfigure.loglog(t[nonZeroDh], dh)

        self.subfigure.set_xlabel("t")
        self.subfigure.set_ylabel("W")

        self.subfigure3.plot(dh, hw[nonZeroDh], 'kx', markersize=0.5)
        self.subfigure3.set_xlabel("RMS(h)")
        self.subfigure3.set_ylabel("hw - m(h)")

        nonzeroEq = where(eqC != 0)
        print nonzeroEq

        if not nonzeroEq:
            return

        self.subfigure4.plot(eqC[nonzeroEq])
        self.subfigure4.set_xlabel("cycle")
        self.subfigure4.set_ylabel("cEq")

        nonzero = where(mc != 0)
        if not nonzero:
            return

        mc = mc[nonzero]

        self.subfigure4.plot(mc, label="mean")
        self.subfigure4.legend()


class logmtvec(DCVizPlotter):

    nametag = "logmtvec*"

    armaBin = True

    def plot(self, data):

        self.subfigure.plot(data.data)


class fit_scale(DCVizPlotter):

    nametag = "scale_fits\.arma";

    armaBin = True

    def plot(self, data):

        im = self.subfigure.imshow(data.data)
        self.subfigure.set_xlabel("power")
        self.subfigure.set_ylabel("scale")
        self.figure.colorbar(im)



class ebs_s(DCVizPlotter):

    nametag = "eb_s_(.*)\.arma"

    isFamilyMember = True

    armaBin = True

    hugifyFonts = True

    def plot(self, data):

        for name, eb_s in zip(self.familyFileNames, data):

            eb, s = eb_s
            eb = eb[:]
            s = s[:]

            self.subfigure.plot(eb, s, self.colors.next() + self.markers.next(), fillstyle="none", label=re.findall(self.nametag, name)[0])
        # self.subfigure.plot(eb, exp(-0.5*eb))
        # self.subfigure.plot(eb, exp(exp(-9/4.*eb+1) - eb/2 + 1./10))

        self.subfigure.set_xbound(0)
        self.subfigure.set_xlabel("Eb/kT")
        self.subfigure.set_ylabel(r"$\langle s\rangle/L$")
        self.subfigure.legend()


class SOSanalyze(DCVizPlotter):

    nametag = "analyze_(.+)\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"figure" : "subfigure", "figure_2" : "mean_figure"}

    def stuff(self, data, dirname, label):
        print dirname

        C =  data[self.get_family_index_from_name("analyze_%s_C_values.npy" % dirname)].data
        s0 =  data[self.get_family_index_from_name("analyze_%s_s0_values.npy" % dirname)].data
        r0 =  data[self.get_family_index_from_name("analyze_%s_r0_values.npy" % dirname)].data

        r0_mean = C.sum(axis=0)/len(s0)

        self.mean_figure.plot(r0, r0_mean, "k^", label=label, linewidth=3,
                                fillstyle='none',
                                markeredgewidth=1.5,
                                markersize=5)

    def plot(self, data):

        dirname = re.findall("analyze\_(.+)\_.+\_values\.npy", self.familyHead)[0]

        print dirname

        C =  data[self.get_family_index_from_name("analyze_%s_C_values.npy" % dirname)].data
        s0 =  data[self.get_family_index_from_name("analyze_%s_s0_values.npy" % dirname)].data
        r0 =  data[self.get_family_index_from_name("analyze_%s_r0_values.npy" % dirname)].data

        s0_cut = where(s0 >= 0)
        r0_cut = where(r0 > 0.4)


        s0_idx, r0_idx = numpy.meshgrid(s0_cut, r0_cut)


        ax = Axes3D(self.figure)


        X, Y = numpy.meshgrid(s0[s0_cut], r0[r0_cut])

        Z = C[s0_idx, r0_idx]

        colormap = cm.hot

        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=colormap, linewidth=0)

        zdir = "x"
        offset = -0.5

        cset = ax.contour(X, Y, Z, zdir=zdir, offset=offset*1.05, color='#008000')

        ax.set_xlim(-0.5, s0[s0_cut].max()*1.1)
        ax.set_zlim(0, Z.max())

        ax.set_xlabel(r"$\sigma_0$")
        ax.set_ylabel(r"$r_0$")
        ax.view_init(30, -65)

        self.subfigure.axes.contourf(X, Y, Z, zdir='z', cmap=colormap)
        self.subfigure.set_xlabel(r"$\sigma_0$")
        self.subfigure.set_ylabel(r"$r_0$", rotation=0)
        # self.subfigure.axes.get_yaxis().get_label().set_fontsize(20)
        # self.subfigure.axes.get_xaxis().get_label().set_fontsize(20)

        r0_mean = C.sum(axis=0)/len(s0)
        r0_mean_fit = r0_mean[r0_cut]

        r0_fit = r0[r0_cut]


        if "noshadow" in dirname:
            label="basic"
            otherlabel="shadowing"

            otherdirname = re.findall("analyze\_(.+)\_.+\_values\.npy", self.familyHead.replace("noshadow", "shadow"))[0]

        else:
            label="shadowing"
            otherlabel="basic"

            otherdirname = re.findall("analyze\_(.+)\_.+\_values\.npy", self.familyHead.replace("shadow", "noshadow"))[0]

        try:
            self.stuff(data, otherdirname, otherlabel)

        except:
            print "error loading: ", otherdirname

        self.mean_figure.plot(r0, r0_mean, 'ks',
                              label=label,
                              linewidth=3,
                              fillstyle='none',
                              markersize=5,
                              markeredgewidth=1.5)

        from scipy.optimize import curve_fit


        f = lambda x, a, b, c: a*x*(1 - exp(-1/x)) + b

        popt, pcov = curve_fit(f, r0_fit, r0_mean_fit, p0=(1., 0., 1))

        a, b, c = popt
        a = 1; b = 0; c = 1

        f = np.vectorize(f)

        self.mean_figure.plot(r0, f(r0, a, b, c), "r-", label="analytical")
        self.mean_figure.legend(loc="lower right")

        self.mean_figure.set_xlabel(r"$r_0$")
        self.mean_figure.set_ylabel(r"$g(r_0)$")
        self.mean_figure.set_xbound(0)
        self.mean_figure.set_ybound(0)

        print a, b, c




class shifts(DCVizPlotter):

    nametag = "shifts\.arma|values\.arma"

    isFamilyMember = True

    armaBin = True

    hugifyFonts = True

    def plot(self, data):

        cut = 5

        # shifts = data[self.get_family_index_from_name("shifts.arma")][0]
        values = data[self.get_family_index_from_name("values.arma")][0]

        print sum(values[cut:])/(len(values) - cut)

        n = avg1 = avg2 = 0
        for v in values[cut:]:
            avg1 += v
            avg2 += v*v
            n+=1
        std = sqrt(1./(n-1)*avg2 - avg1*avg1/(n*(n-1)))
        print std

        # self.subfigure.plot(-shifts,'^', label="shifts")
        self.subfigure.plot(values, "k--s",
                            markersize=7,
                            markeredgewidth=1.5,
                            linewidth=1,
                            label="values",
                            fillstyle="none")

        self.subfigure.plot([0 for i in range(len(values))], 'r--', linewidth=3, label="exact")
        self.subfigure.axes.set_xticks(range(cut + 1))

        self.subfigure.set_xlabel(r"$k$")
        self.subfigure.set_ylabel(r"$\gamma_k$")
        self.subfigure.axes.set_xlim(-0.1, cut)
        self.subfigure.axes.set_ylim(-0.05, 1.1)
        # self.subfigure.legend()

        # zoomfac = len(values)/2
        #
        # if len(shifts) <= zoomfac + 3:
        #     return
        #
        # X = range(zoomfac, len(shifts))
        #
        # self.zoomed.plot(X, shifts[zoomfac:], marker='x', label="shifts")
        # self.zoomed.plot(X, values[zoomfac:], marker='o', label="values")
        # self.zoomed.plot(X, [0 for i in range(len(values) - zoomfac)], 'r--', label="exact")
        # self.zoomed.set_xlabel("n")
        #
        # avg = 0
        # for start in X:
        #     shift_sum = 0
        #     avg2 = 0
        #     for i, value in enumerate(values[start:]):
        #         avg2 += value
        #         shift_sum += 1
        #     avg2 /= shift_sum
        #
        #     print sum(values[start:])/(len(values) - start), avg2
        #
        #     avg += avg2
        #
        #
        # print avg/len(X), values[-1], sum(values)/len(values)



        # self.zoomed.legend()


class s_vs_eb_collapse(DCVizPlotter):

    nametag = "varl_boltzmann_ascii(\d+)\.arma"

    isFamilyMember = True

    def plot(self, family_data):

        for name, data in zip(self.familyFileNames, family_data):

            L = self.getNumberForSort(name)
            ebs, s = data.data

            c = self.colors.next()
            m = self.markers.next()

            self.subfigure.plot(ebs, s, c + m, label="L="+str(L))

        self.subfigure.legend()
        self.subfigure.set_xlabel("alpha")
        self.subfigure.set_ylabel("s/L")

class canonical_stuff(DCVizPlotter):

    nametag = "canonical_(.*)\.arma"

    armaBin = True

    isFamilyMember = True

    hugifyFonts = True

    figMap = {"fig1" : "E_fig",
              "fig2" : "VarE_fig",
              "fig3" : "CV_fig",
              "fig4" : "S_fig",
              "fig5" : "A_fig",
              "fig6" : "logZ_fig"}

    conversion = {"avgS" : "E_fig",
                  "varS" : "VarE_fig",
                  "CV"   : "CV_fig",
                  "S"    : "S_fig",
                  "A"    : "A_fig",
                  "logZ" : "logZ_fig"}

    labels = {"avgS" : r"$\langle E \rangle/E_b$",
              "varS" : r"$\sigma (E) / E_b$",
              "CV"   : r"$C_V/k$",
              "S"    : r"$S/k$",
              "A"    : r"$A/E_b$",
              "logZ" : r"$\log(Z)$"}

    limits = {"avgS" : [[0, None], [0, 10]],
              "varS" : [[0, None], [0, 0.6]],
              "CV"   : [[0, None], [0, None]],
              "S"    : [[0, None], [0, None]],
              "A"    : [[0, None], [-1000, None]],
              "logZ" : [[0, None], [0, None]]}


    def plot(self, family_data):

        alpha_index = self.get_family_index_from_name("canonical_alphas.arma")
        alphas = family_data[alpha_index]

        for i, data in enumerate(family_data):

            if i == alpha_index:
                continue

            which = self.get_nametag_match(self.familyFileNames[i])

            fig = eval("self.%s" % self.conversion[which])

            fig.plot(alphas.data, data.data, "r-", linewidth=3)
            fig.axes.set_xticks(range(int(round(alphas.data.max())) + 1))
            fig.set_xlabel(r"$\alpha$")
            fig.set_ylabel(self.labels[which])

            if self.limits[which]:

                if self.limits[which][0]:
                    xmin, xmax = self.limits[which][0]

                    if not xmin:
                        xmin = fig.axes.get_xlim()[0]
                    if not xmax:
                        xmax = fig.axes.get_xlim()[1]

                    fig.axes.set_xlim([xmin, xmax])

                if self.limits[which]:
                    ymin, ymax = self.limits[which][1]

                    if not ymin:
                        ymin = fig.axes.get_ylim()[0]
                    if not ymax:
                        ymax = fig.axes.get_ylim()[1]

                    fig.axes.set_ylim([ymin, ymax])


class boltzmann_ascii(DCVizPlotter):

    nametag = "boltzmann_ascii(.*)\.arma"

    transpose = True

    hugifyFonts = True

    def plot(self, data):

        alpha, E = data
        logE = log(E)

        self.subfigure.plot(alpha, logE, "r-", linewidth=3)
        self.subfigure.axes.set_xticks(range(int(round(alpha.max())) + 1))
        self.subfigure.set_xbound(0)
        self.subfigure.set_xlabel(r"$\alpha$")
        self.subfigure.set_ylabel(r"$\log\left(\langle s_\uparrow \rangle / L\right)$")

        self.subfigure.plot([0, 1], [0, 0], "k--")
        self.subfigure.plot([1, 1,], [self.subfigure.axes.get_ylim()[0], 0], "k--")

class s_vs_eb_blocked(DCVizPlotter):

    nametag = ".*_blocking\.dat"
    skipCols = 1

    transpose = True

    hugifyFonts = True

    def plot(self, data):

        s, err = data
        err *= 1

        alphas = []
        for name in self.skippedCols[0]:
            alpha = re.findall("alpha_(.+?)_mu", name)[0]
            L = re.findall("L_(.+?)_alpha", name)[0]
            alphas.append(float(alpha))

        analytical_path = os.path.join("/tmp", "boltzmann_ascii%s.arma" % L)
        if os.path.exists(analytical_path):
            analytical = numpy.loadtxt(analytical_path)
            self.subfigure.plot(analytical[:, 0], analytical[:, 1], 'r-',
                                linewidth=3,
                                label="$\mathrm{Analytical}$",
                                fillstyle='none',
                                markersize=5)

        # self.subfigure.axes.errorbar(alphas, s,
        #                              yerr=err,
        #                              fmt='s',
        #                              fillstyle='none',
        #                              label="$\mathrm{KMC}$",
        #                              markersize=7,
        #                              markeredgewidth=1.5,
        #                              linewidth=1,
        #                              color="black")

        self.subfigure.plot(alphas, s, 'ks',
                             fillstyle='none',
                             label="$\mathrm{KMC}$",
                             markersize=7,
                             markeredgewidth=1.5,
                             linewidth=1)

        self.subfigure.legend()

        self.subfigure.set_xlabel(r"$\alpha$")
        self.subfigure.set_ylabel(r"$\langle E \rangle / L$")
        self.subfigure.set_xlim(0, max(alphas))
        self.subfigure.set_ylim(0, max(s)*1.2)

class Quasi2D_slopes_and_stuff(DCVizPlotter):

    nametag = "linearplots\_(.*)\_?\d*\.npy"

    numpyBin = True

    isFamilyMember = True

    figMap = {"fig" : "gammaslopes", "fig2" : "E0slopes"}

    hugifyFonts = True

    def n_runs(self):

        n_max = -1

        for name in self.familyFileNames:
            if not "alpha" in name:
                continue

            n = int(re.findall("linearplots_alpha_(\d+)\.npy", name)[0])

            if n > n_max:
                n_max = n

        return n_max + 1


    def plot(self, data):

        meta_data_file = os.path.join(self.familyHome, "linearplots_metadata.dat")

        with open(meta_data_file, 'r') as meta_data:
            print meta_data.read()

        N = 3
        shapes = ["s", "^", "v"]

        n_runs = self.n_runs()

        E0s = data[self.get_family_index_from_name("linearplots_E0.npy")].data
        slopes = data[self.get_family_index_from_name("linearplots_slopes.npy")].data

        n_plots = 0
        errors = []
        for n in range(n_runs):

            alpha_name = "linearplots_alpha_%d.npy" % n
            mu_name = "linearplots_muEqs_%d.npy" % n

            alphas = data[self.get_family_index_from_name(alpha_name)].data
            mus = data[self.get_family_index_from_name(mu_name)].data

            a, b, c, d, err = linregress(alphas, mus)

            errors.append(err)

            if n%(n_runs/(N-1)) != 0:
                continue
            print n, n_plots

            mu_error_name = "linearplots_muEqErrors_%d.npy" % n

            mu_errors = data[self.get_family_index_from_name(mu_error_name)].data

            self.gammaslopes.errorbar(alphas, mus,
                                      yerr=mu_errors,
                                      fmt=shapes[n_plots],
                                      fillstyle='none',
                                      label=r"$|E_0|/L=%1.2f$" % (E0s[n]),
                                      markersize=7,
                                      markeredgewidth=1.5,
                                      linewidth=1,
                                      color="black")
            self.gammaslopes.plot([0, alphas.max()], [0, slopes[n]*alphas.max()], "r--")

            n_plots += 1

        self.gammaslopes.set_xbound(0)
        self.gammaslopes.set_ybound(0)
        self.gammaslopes.legend(loc="upper left")

        self.gammaslopes.set_xlabel(r"$\alpha$")
        self.gammaslopes.set_ylabel(r"$\gamma_\mathrm{eq}$")

        self.E0slopes.errorbar(E0s, slopes,
                               yerr=errors,
                               fmt="s",
                               fillstyle='none',
                               markersize=7,
                               markeredgewidth=1.5,
                               linewidth=1,
                               color="black")

        sslope, a, b, c, d = linregress(E0s, slopes)

        self.E0slopes.plot([0, E0s.max()], [0, sslope*E0s.max()], 'r--')
        self.E0slopes.set_xbound(0)
        self.E0slopes.set_ylim(0, sslope*E0s.max()*1.05)

        self.E0slopes.set_xlabel(r"$|E_0|/L$")
        self.E0slopes.set_ylabel(r"$\frac{\partial}{\partial \alpha}\gamma_\mathrm{eq}$")

        print sslope, d

class SOS_pressure_sizes(DCVizPlotter):

    nametag = "pressure_plots_.*\.npy"

    numpyBin = True

    isFamilyMember = True

    hugifyFonts = True

    def plot(self, data):

        N = 3
        shapes = ["s", "^", "v"]


        E0_array = data[self.get_family_index_from_name("pressure_plots_E0.npy")].data
        alphas = data[self.get_family_index_from_name("pressure_plots_alphas.npy")].data
        mean_s = data[self.get_family_index_from_name("pressure_plots_mean_s.npy")].data
        var_s = data[self.get_family_index_from_name("pressure_plots_var_s.npy")].data

        print E0_array.shape, alphas.shape, mean_s.shape, var_s.shape

        analytical_path = os.path.join("/tmp", "boltzmann_ascii256.arma")
        if os.path.exists(analytical_path):
            analytical = numpy.loadtxt(analytical_path)
            self.subfigure.plot(analytical[:, 0], analytical[:, 1], 'r-',
                                linewidth=3,
                                label="$\mathrm{Analytical}$",
                                fillstyle='none',
                                markersize=5)

        for i, E0_value in enumerate(E0_array):


            if i%(len(E0_array)/(N-1)) != 0:
                continue

            alpha_array = alphas[i, :]
            mean_s_array = mean_s[i, :]
            var_s_array = var_s[i, :]


            self.subfigure.plot(alpha_array, mean_s_array, 'ks',
                                 fillstyle='none',
                                 label="$E_0 = %g$" % E0_value,
                                 markersize=7,
                                 markeredgewidth=1.5,
                                 linewidth=1)

        self.subfigure.set_xbound(alpha_array.min()*0.9)
        self.subfigure.legend()




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
    
    nametag = "blocking_.*?\.dat"
    figMap = {"Fig": ["blockFig"]}
    
    nameMap = {"0": r"$\alpha$", "1": r"$\beta$", "": ""}
    
    transpose = True    
    
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
    
    nametag = "DMC_out.*\.dat"
    figMap = {"Fig": ["N_plot"], "Fig2": ["E_plot"]}
    dt = 0.001
        
    transpose = True
    
    def plot(self, data):
        print data.data.shape
        
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
    
    cut = 4
    yMax = None
    silent = False
    
    color = ['#008000', "0.75", "k"]
    style = ['-', '-.', '--']

    pureOnly = False

    def plot(self, data):
        
        
        r2 = False
        if "Atoms" in self.familyFileNames[0]:
            r2 = True

        if len(self.familyFileNames) == 1:
            maxCut, max_edge = self.plotsingle(data, 0)
        else:

            root_id = re.findall("radial\_out\_(.*)\_edge", self.originalFilename)[0].replace("dmc", "qmc").replace("vmc", "qmc")
            
            k = None        
            o = None
            for i, name in enumerate(self.familyFileNames):

                if name == self.originalFilename:
                    o = i
                    if k is not None:
                        break
                    else:
                        continue
            
                _id = re.findall("radial\_out\_(.*)\_edge", name)[0].replace("dmc", "qmc").replace("vmc", "qmc")
                
                if (_id == root_id and k is None):
                    k = i
                    if o is not None:
                        break
                
    
            if k is None:
                if not self.silent: print "No matching ids found."
                
                maxCut, max_edge = self.plotsingle(data, o)                
                
            else:
                if not self.silent: print self.familyFileNames[o], self.familyFileNames[k]
                
                maxCut, max_edge = self.plotdouble(data, o, k)
               
            
        self.radialFig.legend()    
        self.radialFig.axes.set_xlim(maxCut, max_edge)
        self.radialFig.set_xlabel('r')
        
        if r2:
            self.radialFig.set_ylabel(r'$r^2\rho(r)$')
        else:
            self.radialFig.set_ylabel(r'$\rho(r)$')

        if self.yMax is not None:
            self.radialFig.set_ylim(0, self.yMax)
        self.radialFig.axes.set_ybound(0)
        
        self.radialFig.axes.get_yaxis().get_label().set_fontsize(30)
        self.radialFig.axes.get_xaxis().get_label().set_fontsize(30)

    def getEdge(self, n):
        return float(re.findall("edge(\d+\.?\d*)\.arma", self.familyFileNames[n])[0])
          
    def plotdouble(self, data, o, k):
        
        cuts = []
        edges = []
        
        if not self.pureOnly:
            cut0, edge0 = self.plotsingle(data, o)
            cut1, edge1 = self.plotsingle(data, k, 1)
            
            cuts.append(cut0)
            cuts.append(cut1)
            
            edges.append(edge0)
            edges.append(edge1)
            
        if "vmc" in self.familyFileNames[o] and "dmc" in self.familyFileNames[k]:
            vmcDist = data[o].data
            dmcDist = data[k].data
        elif "vmc" in self.familyFileNames[k] and "dmc" in self.familyFileNames[o]:
            vmcDist = data[k].data
            dmcDist = data[o].data
        else:
            raise RuntimeError("Matching distribution sets are not vmc dmc pairs.")
      
        dist = 2*dmcDist - vmcDist
        
        edge = self.getEdge(o)
        
        r = linspace(0, edge, data[o].n)
        
        self.radialFig.plot(r, dist, self.style[2], color=self.color[2], label="Pure")
        
        cuts.append(r[self.cut])
        edges.append(edge)
                
        
        return max(cuts), max(edges)        
        
    def plotsingle(self, data, n, i = 0):
        edge = self.getEdge(n)
        
        r = linspace(0, edge, data[n].n)

        if "vmc" in self.familyFileNames[n]:
            label = "VMC"
        else:
            label = "DMC"
            
        self.radialFig.plot(r, data[n].data, self.style[i], color=self.color[i], label=label)


        if not self.silent: print self.familyFileNames[n]
        if not self.silent: print "n_p= ", data[n].data.sum()*edge/(data[n].n-1), "?"
            
        
        return r[self.cut], edge
        
        

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
    
        if len(data) == 1:
            if not silent: print "length 1 data"
            
            edge = float(re.findall("_edge(.+?)\.arma", self.familyFileNames[0])[0])
            dist = data[0].data
            
        else:

            root_id = re.findall("dist\_out\_(.*)\_edge", self.originalFilename)[0].replace("dmc", "qmc").replace("vmc", "qmc")
            
            k = None        
            o = None
            for i, name in enumerate(self.familyFileNames):

                if name == self.originalFilename:
                    o = i
                    if k is not None:
                        break
                    else:
                        continue
            
                _id = re.findall("dist\_out\_(.*)\_edge", name)[0].replace("dmc", "qmc").replace("vmc", "qmc")
                
                if (_id == root_id and k is None):
                    k = i
                    if o is not None:
                        break
            

            edge = float(re.findall("_edge(.+?)\.arma", self.familyFileNames[o])[0])            
            
            if k is not None:
                if not silent: print self.familyFileNames[o], self.familyFileNames[k]
                edge_2 = float(re.findall("_edge(.+?)\.arma", self.familyFileNames[k])[0])
            
                if not silent:
                    if edge != edge_2:
                        print "Bin edges does not match. %s != %s" % (edge, edge_2)
    
                if "vmc" in self.familyFileNames[o] and "dmc" in self.familyFileNames[k]:
                    vmcDist = data[o].data
                    dmcDist = data[k].data
                elif "vmc" in self.familyFileNames[k] and "dmc" in self.familyFileNames[o]:
                    vmcDist = data[k].data
                    dmcDist = data[o].data
                else:
                    raise RuntimeError("Matching distribution sets are not vmc dmc pairs.")
                
                if self.dmcOnly:
                    try:
                        dist = dmcDist
                    except: 
                        if not silent: print "\n\nWarning: No DMC data found. Attempting to load VMC data.\n\n"
                       
                        dist = vmcDist
                       
                elif self.vmcOnly:
                    dist = vmcDist
                else:
                    dist = 2*dmcDist - vmcDist
                    print "pure success!", dmcDist.sum()*(edge/100)**2, vmcDist.sum()*(edge/100)**2 
            
            else:
                if not silent: print "length 1 data"
                
                dist = data[o].data
        
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
        
        try:
            from scipy.stats import linregress as l2
        except:
            pass
            
              
        
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

class HeatCap(DCVizPlotter):
    
    nametag = "CT\.arma"
    armaBin = True    
    
    def plot(self, data):
        
        C, T = data        
        
        self.subfigure.plot(T, C)
        self.subfigure.set_ybound(0)
        self.subfigure.set_ylabel("C(T)")
        self.subfigure.set_xlabel("T")
        
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

            try:
                from scipy.stats import linregress as l2
            except:
                pass

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
    
    nametag = "ASGD_out.*\.dat"
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
