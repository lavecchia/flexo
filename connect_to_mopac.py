import objects
import commonfunctions
from commonparameters import *

#read output of MOPAC
def mopac_outputread(namefile): 
    heatformation = None
    atomlst=[]
    mopacfile = open(namefile,"r")
    lines = mopacfile.readlines()
    mopacfile.close()
    indexcoord=lines.index("                             CARTESIAN COORDINATES\n")
    indexcoord += 2
    while lines[indexcoord]!="\n":
        linediv = lines[indexcoord].split()
        symbol = linediv[1]
        number = int(linediv[0])
        xcor = float(linediv[2])
        ycor = float(linediv[3])
        zcor = float(linediv[4])
        atom = objects.myAtom(symbol,number,xcor,ycor,zcor)
        atomlst.append(atom)
        indexcoord += 1
    
    for line in lines:
        if KEYFHOF in line:
            linediv = line.split()
            #extract FINAL HEAT OF FORMATION in KCAL/MOL
            heatformation = float(linediv[5])
        elif KEYCA in line:
            linediv = line.split()
            #extract COSMO AREA in SQUARE ANGSTROMS
            cosmoarea = float(linediv[3])
    
    if heatformation:
        return heatformation,cosmoarea,atomlst
    else:
        # If there is a problem, perhap there is a clash
        print "ERROR: A mistake produced to try read the file " + namefile
        return None, None, None
   

#make MOPAC input file from a molecule object
def mopac_makeinput(outputfilename, conformer, keywords=""):
    outputfile = open(outputfilename,"w")
    outputfile.write("%s \n\n\n"%(keywords))
    for atom in conformer.atomlst:
        # check if the atom must be freezed
        if atom.freeze == True:
            tabfreeze = "0"
        else:
            tabfreeze = "1"
        outputfile.write("%s  %8.5f %s %8.5f %s %8.5f %s\n"%(atom.symbol,atom.xcor, tabfreeze, atom.ycor, tabfreeze, atom.zcor, tabfreeze))
    outputfile.close()
    return 0
    

def mopac_runjobs(conformerlst, keywords="",opt=False):
    if opt==True:
        if "1SCF" in keywords:
            keywords=keywords.replace("1SCF","OPT")
        elif "OPT" in keywords:
            pass
        else:
            keywords+= " OPT"
    
    conformercount = 0
    mopacjoblst = []
    for conformer in conformerlst:
        conformercount += 1
        inputfilename = "%i.mop"%(conformercount)
        mopac_makeinput(inputfilename, conformer, keywords)
        mopacjoblst.append(inputfilename)
    commands = []
    for inputfilename in mopacjoblst:
        commands.append([MOPACPATH, inputfilename])
    commonfunctions.exec_commands(commands)        
             
    for conformer,inputfilename in map(None,conformerlst,mopacjoblst):
        hofcosmo, areacosmo, atomlst = mopac_outputread(inputfilename.replace(".mop",".out")) #read MOPAC output from .out
        conformer.energy = hofcosmo
        conformer.sasa = areacosmo
        if opt==True:
            conformer.atomlst = atomlst
        #~ cosmoparamlist = cosmoout_read(cosmoinput.replace(".mop",".cos")) #read cosmo param from .cos
        #~ electronlist = electron_read(cosmoinput.replace(".mop",".out")) #read number of electron by atom from MOPAC output 
        
    return conformerlst
        
     
