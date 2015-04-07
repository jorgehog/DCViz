# -*- coding: utf-8 -*-

import os, sys, inspect, re, threading, time, imp
from os.path import join as pjoin

try:
    import matplotlib
except:
    print "\n"
    print "You need matplotlib in order to run this library!"
    print "sudo apt-get install python-matplotlib"
    print "\n"
    sys.exit(1)

#Adding the srcDir to the local pythonpath in order to avoid global pythonpath sets
srcDir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

def terminalTracker(head, msg):
    print "[%s] %s" % (head.center(10), msg)

def is_DCViz_file(file):
    return file.endswith(".py")

def already_loaded(class_files, file):
    for class_file in class_files:
        if os.path.samefile(class_file, file):
            return True
    return False

def load_class_files():

    try:
        paths = os.environ['DCVIZ_PATH'].split(":")
    except KeyError:
        paths = []

    if srcDir not in paths:
        paths.append(srcDir)

    class_files = []

    for path in paths:
        if os.path.isdir(path):

            for file in os.listdir(path):
                filepath = os.path.join(path, file)
                if os.path.isfile(filepath):
                    if is_DCViz_file(file):
                        if not already_loaded(class_files, filepath):
                            class_files.append(filepath)

        elif os.path.isfile(path):
            if is_DCViz_file(path):
                if not already_loaded(class_files, filepath):
                    class_files.append(path)

        else:
            pass
            # raise RuntimeWarning("DCViz path %s is not eligable as a class file location." % path)

    return class_files

def autodetectModes():

    classfiles = load_class_files()

    if not classfiles:
        terminalTracker("Warning", "No class files found.")

    all_modes = []
    all_modes_names = []
    all_modules = []
    for classfile_name in classfiles:

        classfile = open(classfile_name, 'r')
        raw = classfile.read()
        classfile.close()

        uniqueModesNames = re.findall('^class (\w+)\(DCVizPlotter\):', raw, re.MULTILINE)

        if not uniqueModesNames:
            continue

        terminalTracker("Detector".center(10), "Parsing class file %s" % classfile_name)

        path, name = os.path.split(classfile_name)
        mod_name = name.strip(".py")
        mod = imp.load_source(mod_name, classfile_name)

        uniqueModes = [eval("mod.%s" % subclass) for subclass in uniqueModesNames]


        terminalTracker("Detector".center(10), "Found subclasses %s" \
                                % str(uniqueModesNames).strip("]").strip("["))
        for mode in uniqueModes:

            instance = mode()
            try:
                instance.nametag
            except:
                terminalTracker("Warning", "Subclass %s has no attribute 'nametag' (output filename identifier)." % \
                                  uniqueModesNames[uniqueModes.index(mode)])

        all_modes += uniqueModes
        all_modes_names += uniqueModesNames
        all_modules.append([mod_name, classfile_name])

    return all_modes, all_modes_names, all_modules
    
def matchMode(modes, path, noWarnings=False, silent=False):
    root, name = os.path.split(path)    
    
    if modes is None:
        modes = autodetectModes(path)[0]
    
    matchedMode = None
    for mode in modes:
        if re.findall(mode.nametag, name):
            matchedMode = mode
            break
    
    if matchedMode is None:
        if not noWarnings:
            terminalTracker("Warning", "Found no matching nametags for specified filename %s (%s)" % (name, path))
        return
    
    if not silent: terminalTracker("DCViz", "Matched [%s] with [%s]" % (name, "".join(str(matchedMode).split(".")[1:])))
    return matchedMode

def getInstance(path, dynamic=False, toFile=False, threaded=False):
    
    if dynamic:
        
        N = 0        
        while not os.path.exists(path) and N < 10:
                time.sleep(0.1)
                N += 1
                
    if not os.path.exists(path):
        terminalTracker("Warning", "No such file: " + path)
        return
        
    modes = autodetectModes()[0]
    matchedMode = matchMode(modes, path)
    
    if not matchedMode:
        return
        
    return matchedMode(path, dynamic=dynamic, toFile=toFile, threaded=threaded)

def main(path, app_argv, dynamic, delay = None, toFile=False, fs=None, silent=False, makeGif=False):

    if fs:
        pylab.rcParams['figure.figsize'] = fs

    instance = getInstance(path, dynamic=dynamic, toFile=toFile)

    if not instance:
        return 
        
    if delay is not None and type(delay) in [int, float]:
        instance.delay = delay
        instance.makeGif = makeGif

    instance.mainloop(argv=app_argv)

def mainToFile(path, app_argv):

    name_set = False
    name = ""

    if not os.path.isdir(path):

        if os.path.isfile(path):
            path, name = os.path.split(path)
            name_set = True

            if not path:
                path = os.getcwd()

        else:
            raise Exception("Supplied path must be to a directory or file.")
    
    modes = autodetectModes()[0]
    print path
    for root, dirs, files in os.walk(path):
    
        init = True    
    
        for outfile in sorted(files):

            if name_set:
                if outfile != name:
                    continue

            matchedMode = matchMode(modes, pjoin(root, outfile), noWarnings=True)

            if matchedMode is None:
                continue
            
            if matchedMode.isFamilyMember:
                if init:
                    thisTag = matchedMode.nametag
                    init = False
                else:
                    if re.findall(thisTag, outfile):
                        terminalTracker("DCViz", "Family portait already taken, skipping...")                        
                        continue
                    else:
                        thisTag = matchedMode.nametag
            
            plotTool = matchedMode(pjoin(root, outfile), toFile=True)

            try:
                plotTool.mainloop(argv=app_argv)
            except:
                terminalTracker("DCViz", "Unable to save %s to file: Crash!" % outfile)


class DCVizThread(threading.Thread):
    
    def __init__(self, filepath, argv, dynamic = False, toFile = False, delay = None):
        
        self.app = getInstance(filepath, dynamic, toFile, threaded=True)
        self.argv = argv

        if delay is not None:
            
            if type(delay) not in [float, int]:
                terminalTracker("DCViz", "Delay must be given as a numeric value")                

            self.app.delay = delay
        
        super(DCVizThread, self).__init__()

    def run(self):
        self.app.mainloop(argv=self.argv)
        
    def stop(self):
        
        if not self.app.stopped:
            self.app.stopped = True
            terminalTracker("DCViz", "Stopped thread.")


if __name__ == "__main__":
    dynamic = False
    toFile = False
    makeGif = False
    delay = None

    dcviz_argv = []
    path = None
    for k, arg in enumerate(sys.argv[1:]):
        if os.path.exists(arg) or os.path.exists(os.path.join(os.getcwd(), arg)):
            dcviz_argv = sys.argv[1:k+1]
            app_argv = sys.argv[k+2:]

            path = arg
            break

    if not path:
        raise RuntimeError("Path error.")


    if "-d" in dcviz_argv:
        dynamic = True

        idx = dcviz_argv.index("-d")

        if not idx == len(dcviz_argv):
            obj = eval(dcviz_argv[idx + 1])
            if type(obj) in [int, float]:
                delay = obj
                dcviz_argv.pop(idx+1)

        dcviz_argv.remove("-d")

    if "-g" in dcviz_argv:
        if not dynamic:
            print "Making GIF requires dynamic mode."
            sys.exit(1)

        makeGif=True
        dcviz_argv.remove("-g")

    elif "-f" in dcviz_argv:

        if makeGif:
            print "Either make a gif or dump figures."
            sys.exit(1)

        toFile = True
        dcviz_argv.remove("-f")

    figSize = None
    if "-s" in dcviz_argv:
        idx = dcviz_argv.index("-s")

        if not idx + 1 == len(dcviz_argv) and not idx + 2 == len(dcviz_argv):
            figSize = [eval(dcviz_argv[idx + 1]), eval(dcviz_argv[idx + 2])]
            for obj in figSize:
                if type(obj) in [int, float] and figSize is not None:
                    dcviz_argv.pop(idx+ 1)
                else:
                    figSize = None

        dcviz_argv.remove("-s")


    if len(dcviz_argv) != 0:
        print "Unrecognized command line arguments:"
        for arg in dcviz_argv:
            print arg
        sys.exit(1)

    if toFile:
        mainToFile(path, app_argv)
    else:
        main(path, app_argv, dynamic, delay=delay, fs=figSize, makeGif=makeGif)