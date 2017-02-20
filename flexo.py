# FLEXO Conformational search
# INPUT: 
# mol2 structure

#import guppy
#h=guppy.hpy()
#print h.heap()


import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

from commonfunctions import * # common function
from commonparameters import * # common parameters
import numpy as np
import scipy.optimize
import subprocess
from itertools import count
import ConfigParser # import setting file
from optparse import OptionParser
import sys, os, time
from subprocess import Popen, list2cmdline
import glob
import time #check time of calculation
import datetime
import shutil
import math
import ga_interface
from mainfunctions import *

pymol.finish_launching()

###########
# PARSE
###########
parser = OptionParser()
parser.add_option("-o", "--outfile", dest="outfilename", default=DEFREPORFILE,
                  help="name of file to write REPORT", metavar="REPORT")
parser.add_option("-c", "--config", dest="configfile", default=CONFIG_FILENAME,
                  help="FILE to setting the program", metavar="CONFIG")
parser.add_option("-i", "--infile", dest="infile", 
                  help="Molecular structure in .mol2", metavar="INFILE")

(options, args) = parser.parse_args()

outfilename = options.outfilename
infilename = options.infile


cfg = ConfigParser.ConfigParser() 
if not cfg.read([options.configfile]):  
    print "Running without settings file"  


###########
# VARS
###########
varlist = [] #store name of vars to then print
mcmarklist = [] #store the last mcmark to control the temperature in MC run
#limits
if cfg.has_option("limits", "dihedral_minvalue"):  
    dihedral_minvalue = float(cfg.get("limits", "dihedral_minvalue"))
else:  
    dihedral_minvalue = -179.9
varlist.append("dihedral_minvalue")

if cfg.has_option("limits", "dihedral_maxvalue"):  
    dihedral_maxvalue = float(cfg.get("limits", "dihedral_maxvalue"))
else:  
    dihedral_maxvalue = 180.0
varlist.append("dihedral_maxvalue")


#method
if cfg.has_option("method", "forcefield"):  
    forcefield =cfg.get("method", "forcefield")
else:  
    forcefield = "mopac"
varlist.append("forcefield")

if cfg.has_option("method", "extrakeys"):  
    extrakeys =cfg.get("method", "extrakeys").upper()
else:  
    extrakeys = "PM6 PRECISE 1SCF GEO-OK CHARGE=0"
varlist.append("extrakeys")



#system
if cfg.has_option("system", "parameterfile"):
    parameterfile = int(cfg.get("system", "parameterfile"))
else:
    parameterfile = None
varlist.append("parameterfile")


if cfg.has_option("system", "calculationtype"):  
    calculationtype = cfg.get("system", "calculationtype")
else:  
    calculationtype="ga"
varlist.append("calculationtype")

# genetic algorithm parameters
if cfg.has_option("system", "numbermembers"):
    numbermembers = int(cfg.get("system", "numbermembers"))
else:
    numbermembers = 10
varlist.append("numbermembers")
    
if cfg.has_option("system", "maxgen"):
    maxgen = int(cfg.get("system", "maxgen"))
else:
    maxgen = 10
varlist.append("maxgen")

if cfg.has_option("system", "bestconfnum"):
    try:
        bestconfnum = int(cfg.get("system", "bestconfnum"))
    except:
        bestconfnum = cfg.get("system", "bestconfnum")
    
else:
    bestconfnum = BESTCONFNUM
varlist.append("bestconfnum")



###########
# PROGRAM
###########

# NOTATION SUFFIX: 
#         0: initial value of the structure
#      best: the lowest energy conformers
#   current: store current values of run
print forcefield
#~ mask = "3972:4014"
mask = None

finalopt=True

molecule=myMolecule(infilename)

# PARAMETER TO MODIFY:
# List of parameter objects
if parameterfile==None:  
    dihedrallst=identify_dihedral(molecule,mask)

paramlst = []
for dihedral in dihedrallst:
    #~ print "checking dihedral:"
    #~ print dihedral
    atom1=dihedral[0]
    atom2=dihedral[1]
    atom3=dihedral[2]
    atom4=dihedral[3]
    value=get_dihedral(molecule.filename,atom1,atom2,atom3,atom4)
    paramtype="D@"
    paramlst.append(Parameter(paramtype,dihedral,value,dihedral_minvalue,dihedral_maxvalue))
molecule.initialparamlst=paramlst
if len(paramlst) == 0:
    print "There is not parameters to modify"
    exit()
else:
    print "Parameters to modify: %i"%(len(paramlst))


#initial time
start_time = time.time()

#initialize report file
#~ outfile = open(outfilename, "w",0)
outfile = open(outfilename, "w")

#head of report file
outfile.write(str(datetime.datetime.now())+"\n")

#write settings in report file
for name in varlist:
    try:
        outfile.write(name + "\t\t=\t" + str(eval(name)) + "\n")
    except:
        pass

bestconformers = bestConformerList(bestconfnum)

if calculationtype == 'ga': 
    print "Genetic Algorithm start"
    conformerlst = [] #list that store members or individues
    
    #building parents
    #~ while len(conformerlst)<numbermembers:
        #~ # generate new values of parameters with a distribution probability from current values
        #~ paramtestlst = generate_values(molecule.initialparamlst)
        #~ 
        #~ 
        #~ if check_clashes(newconformer) == False:
            #~ conformerlst.append(molecule.params_to_conformer(paramtestlst))
        #~ else:
            #~ print "A clash was detected"
            #~ del newconformer
            #~ del paramtestlst
   

    #~ conformerlst = calc_energy(forcefield, conformerlst, extrakeys)

    #generations
    nclash=0 
    for ncycle in range(1,maxgen+1):
        startcycle_time = time.time()
     
        #generate a new conformer if one had clashes and was eliminated 
        while len(conformerlst)<numbermembers:
			#use the initial conformer to generate a new member
            newparamlst = generate_values(molecule.initialparamlst)
            newconformer = molecule.params_to_conformer(newparamlst)
            
            #~ del newparamlst #memory problem
            
            if check_clashes(newconformer) == False:
                conformerlst.append(newconformer)
            else:
                #~ print "A clash was detected"
                del newconformer
                del newparamlst
                nclash+=1
                if nclash==100:
                    print "A clash was detected"
                    nclash=0

        #update energy
        conformerlst = calc_energy(forcefield, conformerlst, extrakeys)
        energylst = []
        for conformer in conformerlst:
            energylst.append(conformer.energy)
        print "CYCLE: %i AVERAGE ENERGY: %.3f TIME: %s s"%(ncycle,np.mean(energylst),str(time.time()-startcycle_time))
        bestconformerlst = bestconformers.update_list(conformerlst)
        #first remove xyz output files because bestconformer number change every cycle     
        os.system("rm %s_out*.xyz >/dev/null 2>&1"%(molecule.filename.replace(".mol2","")), )
            
        for i in range(0,len(bestconformerlst)):
            conformer_to_xyz(bestconformerlst[i], molecule.filename.replace(".mol2","_out%i.xyz"%(i)))
       
        genetic = ga_interface.SetsGeneration(conformerlst,molecule)
        conformerlst = genetic.next()
        
    if finalopt==True:
        bestconformeroptlst = calc_energy(forcefield,bestconformerlst,extrakeys,opt=True)
        bestconformersopt = bestConformerList(bestconfnum)
        bestconformeroptlst = bestconformersopt.update_list(bestconformeroptlst)
        #first remove xyz output files because bestconformer number change every cycle     
        os.system("rm %s_out*.xyz >/dev/null 2>&1"%(molecule.filename.replace(".mol2","")), )
        for i in range(0,len(bestconformeroptlst)):
            conformer_to_xyz(bestconformeroptlst[i], molecule.filename.replace(".mol2","_opt_out%i.xyz"%(i)))
    
    
outfile.write(str(datetime.datetime.now()))
outfile.write(" TOTAL TIME:" + str(time.time() - start_time) + "seconds")

outfile.close() #close report file







