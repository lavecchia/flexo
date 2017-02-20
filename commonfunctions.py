# common functions
from commonparameters import *
import numpy as np
import scipy.optimize
import subprocess
from itertools import count
import sys, os, time
from subprocess import Popen, list2cmdline
import glob
import time #check time of calculation
import datetime
import shutil
import math
import genetic
from scipy import stats
from objects import *
import math

#~ 
#~ import __main__
#~ __main__.pymol_argv = [ 'pymol', '-qcr'] # Quiet and no GUI
#~ pymol.finish_launching()

# pymol environment
#~ moddir='/opt/pymol-svn/modules'
#~ sys.path.insert(0, moddir)
#~ os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')

# pymol launching
#~ import pymol
#~ pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
#~ pymol.finish_launching()



# ===============================================================================================
# READ FILES FUNCTIONS
# ===============================================================================================
                
#read parameters from a file and return dictonary with format [atsymbol: parameter value...]
def parameters_read(namefile):
    paramdic = {}
    try:
        file = open(namefile,"r")
    except:
        print "FILE ERROR: %s has some problems or is not defined. Check it please" %(namefile)
    for line in file:
        print line
        linediv = line.split()
        paramdic[linediv[0]] = float(linediv[1]) 
    file.close()
    return paramdic


# ===============================================================================================
# PRINT AND CONVERSION FUNCTIONS
# ===============================================================================================
    
#transform symbol to atomicnumber
def symbol_to_atomicnumber(symbol):
    try:
        return PERIODICTABLE[symbol]
    except:
        print "ERROR: the symbol " + symbol + " is not setting in function symbol_to_atomicnumber. Please check this"
        exit    
    

#transform atomicnumber to symbol
def atomicnumber_to_symbol(atomicnumber):
    try:
        for symbol, Z in PERIODICTABLE.iteritems():
            if atomicnumber == Z:
                return symbol
    except:
        print "ERROR: the symbol " + str(atomicnumber) + " is not setting in function atomicnumber_to_symbol. Please check this"
        exit()


#print parameters dictionary to text
def print_param(paramdic):
    line=""
    for key, value in sorted(paramdic.iteritems()):
        line += " %-2s % .5f " % (key,value)
    line = line.strip()
    return line
    
    
def print_summary(totalerror, mae, rmse, bias, r2, slope, intercept,ncycle,datadic,nptype):
    summaryfile = open(DEFSUMMARYFILE,"w")
    summaryfile.close()


#generate list of files inside directory with extension ext
def list_files(directory,ext=DEFINPUTEXT):
    filelist = glob.glob(directory + "/*" + ext)
    return sorted(filelist)

#run calculation in serial or parallel
def exec_commands(cmds, cores = DEFCORE):
    ''' Exec commands in parallel in multiple process 
    (as much as we have CPU)
    '''
    if not cmds: return # empty list

    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)

    max_task = cores
    processes = []
    while True:
        while cmds and len(processes) < max_task:
            task = cmds.pop()
            #print list2cmdline(task)
            processes.append(Popen(task,shell=False,
    stdout=subprocess.PIPE,  stderr=subprocess.PIPE))

        for p in processes:
            if done(p):
                if success(p):
                    processes.remove(p)
                else:
                    fail()

        if not processes and not cmds:
            break
        else:
            pass
            #~ time.sleep(1E-13)
            

def calculate_distance(punto1,punto2):
     dist = math.sqrt((punto1[0] - punto2[0])**2 + (punto1[1] - punto2[1])**2 + (punto1[2] - punto2[2])**2)
     return dist
     
     
    
