# -*- coding: utf-8 -*-

import os, sys, inspect, re, threading
from os.path import join as pjoin

#Adding the srcDir to the local pythonpath in order to avoid global pythonpath sets
srcDir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.append(srcDir)
from DCViz_classes import *


def terminalTracker(head, msg):
    print  "[%s] %s" % (head.center(10), msg)

def autodetectModes():
    classfile = open(os.path.join(srcDir, 'DCViz_classes.py'), 'r')
    raw = classfile.read()
    classfile.close()

    uniqueModesNames = re.findall('^class (\w+)\(DCVizPlotter\):', raw, re.MULTILINE)
    uniqueModes = [eval(subclass) for subclass in uniqueModesNames]

  
    terminalTracker("Detector".center(10), "Found subclasses %s" \
                            % str(uniqueModesNames).strip("]").strip("["))


    if not uniqueModesNames:
        terminalTracker("Warning", "No subclass implementations found.")

    for mode in uniqueModes:
        instance = mode()
        try:
            instance.nametag
        except:
            terminalTracker("Warning", "Subclass %s has no attribute 'nametag' (output filename identifier)." % \
                              uniqueModesNames[uniqueModes.index(mode)])
                              
    return uniqueModes
    
def matchMode(modes, path, noWarnings=False, silent=False):
    root, name = os.path.split(path)    
    
    if modes is None:
        modes = autodetectModes(path)
    
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
    
    if not os.path.exists(path):
        terminalTracker("Warning", "No such file: " + path)
        return
        
    modes = autodetectModes()
    matchedMode = matchMode(modes, path)
    
    if not matchedMode:
        return;
        
    return matchedMode(path, dynamic=dynamic, toFile=toFile, threaded=threaded)

def main(path, dynamic, delay = None, toFile=False, silent=False):

    instance = getInstance(path, dynamic=dynamic, toFile=toFile)

    if not instance:
        return 
        
    if delay is not None and type(delay) in [int, float]:
        instance.delay = delay
    
    instance.mainloop()

def mainToFile(path):
    
    if not os.path.isdir(path):
        raise Exception("Supplied path must be to a directory.")
    
    modes = autodetectModes()
    for root, dirs, files in os.walk(path):
    
        init = True    
    
        for outfile in sorted(files):
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
            plotTool.mainloop()
            

class DCVizThread(threading.Thread):
    
    def __init__(self, filepath, dynamic = False, toFile = False, delay = None):
        
        self.app = getInstance(filepath, dynamic, toFile, threaded=True)
        

        if delay is not None:
            
            if type(delay) not in [float, int]:
                terminalTracker("DCViz", "Delay must be given as a numeric value")                

            self.app.delay = delay
        
        super(DCVizThread, self).__init__()
    
    
    def run(self):
        self.app.mainloop()
        
    def stop(self):
        
        if not self.app.stopped:
            self.app.stopped = True
            terminalTracker("DCViz", "Stopped thread.")


if __name__ == "__main__":
    dynamic = False
    toFile = False   
    delay = None
    
    try:
        if "-d" in sys.argv:
            dynamic = True
        
            idx = sys.argv.index("-d")
            
            if not idx + 1 == len(sys.argv):
                obj = eval(sys.argv[idx + 1])
                if type(obj) in [int, float]:
                    delay = obj
                    sys.argv.pop(idx+1)                    
                                
                    
            sys.argv.remove("-d")
            
        elif "-f" in sys.argv:
            toFile = True
            sys.argv.remove("-f")
            
        path = sys.argv[1]
        sys.argv.pop(1)
            
    except:
        print "Please supply a path as cmdline arg"
        sys.exit(1)

    if len(sys.argv) != 1:
        print "Unrecognized command line arguments:"
        for arg in sys.argv[1:]:
            print arg
        sys.exit(1)
  
    if toFile:
        mainToFile(path)
    else:
        main(path, dynamic, delay=delay)