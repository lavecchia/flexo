# ===============================================================================================
# MAIN FUNCTIONS
# ===============================================================================================
   

from commonparameters import *
from connect_to_mopac import *
from connect_to_openbabel import *
import commonfunctions 

import numpy as np
import os
import time
import copy
from pymol import cmd
import pymol

#uncomment to debug
#from guppy import hpy 

def calc_energy(forcefield,conformerlst,keywords,opt=False):
    tocalclst = []
    calculatedlst = []
    
    for conformer in conformerlst:
        #check energy of conformer. If opt option is True, always optimize
        if conformer.energy==None or opt==True:
            tocalclst.append(conformer)
        else:
            calculatedlst.append(conformer)
    if forcefield=="mopac":
        conformerlst = calculatedlst + mopac_runjobs(tocalclst,keywords,opt)
        #uncomment to debug
        #~ memory = hp.heap()
        #~ print memory
        return conformerlst
        
    elif forcefield=="mm":
        conformerlst = calculatedlst + conformer_to_mmenergy(tocalclst, keywords)
        #uncomment to debug
        #~ memory = hp.heap()
        #~ print memory
        return conformerlst
        
    else:
        print "You must define a forcefield"
        exit()
    


#  
#  name: make_gaussmodification
#  Return an aleatory number with a Gaussian distribution obtained from a number. This number is inside a range determined by rangevalue.
#
#  @param value0: Number where the Gaussian distribution is center 
#  @type value0: number
#
#  @param rangevalue: Specify the modification range.
#  @type rangevalue: number
#
#  @return Modified value.
#  @rtype  number
#
def make_gaussmodification(value0,rangevalue):
    #gauss
    #~ return value0 + rangevalue * np.random.normal(0,0.2)
    
    #uniform
    return value0 + rangevalue * np.random.uniform(-1,1)


#  
#  name: modified_values
#  Return a new list with modified values
#
#  @param paramdic: contain name and value of parametes to optimize
#  @type paramdic: dictonary
#
#  @param freeparamlist: name of parameters free of modified
#  @type rangevalue: list
#
#  @return rrange, gammarange, rsolvrange
#  @rtype  number
#
def modified_values(paramlst):
    newparamlst = [] #store new values of parameters
    for param in paramlst:
        rangevalue = (param.maxvalue-param.minvalue)/15.0
        newparam = copy.copy(param) # problem memory
        #~ newparam = param
        newparam.value=make_gaussmodification(param.value,rangevalue)
        newparamlst.append(newparam)
	paramlst = newparamlst
    return paramlst
    

def generate_values(paramlst):
    newparamlst = [] #store new values of parameters
    for param in paramlst:
        newparam = copy.copy(param) # problem memory
        #~ newparam = param
        newparam.value = np.random.uniform(param.minvalue,param.maxvalue)
        newparamlst.append(newparam)
    #~ newparamlst = modified_values(paramlst)
    paramlst = newparamlst
    return paramlst
   

def fit_lineal(x,y):  
    slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x, y)
    return slope, intercept, r_value*r_value
    
            
def identify_rotatable(moleculefilename,mask=None):
    os.system("perl %s %s  >/dev/null 2>&1"%(FINDROTATABLEPATH, moleculefilename))
    rotdatafile = open(moleculefilename.replace(".mol2",".aux"),"r")
    rotatablelst = []
    for line in rotdatafile:
        if "ROT_BOND" in line:
            linediv=line.split()
            atom1 = int(linediv[1])
            atom2 = int(linediv[2])
            if mask!=None:
                minatom = int(mask.split(":")[0])
                maxatom = int(mask.split(":")[1])
                if (minatom<= atom1 <= maxatom) and (minatom<= atom2 <= maxatom):
                    rotatablelst.append([atom1,atom2])
            else:
                rotatablelst.append([atom1,atom2])
    rotdatafile.close()
    return rotatablelst
    

def identify_rotatable2(molecule,mask=None):
    bond1lst = []
    rotatablelst = []
    connectionnumberlst = [0] * len(molecule.initialconformer.atomlst)
    # find terminal atoms
    for atom1, atom2, bondorder in molecule.bondlst:
        connectionnumberlst[atom1-1]+=1
        connectionnumberlst[atom2-1]+=1
    
       
    
    for atom1, atom2, bondorder in molecule.bondlst:
        #bond order is string
        try:
            if (int(bondorder)==1) and (connectionnumberlst[atom1-1]>1) and (connectionnumberlst[atom2-1]>1):
                if molecule.isinringlst[atom1-1]==False or molecule.isinringlst[atom2-1]==False:
                    if mask!=None:
                        minatom = int(mask.split(":")[0])
                        maxatom = int(mask.split(":")[1])
                        if (minatom<= atom1 <= maxatom) and (minatom<= atom2 <= maxatom):
                            
                            rotatablelst.append([atom1,atom2])
                    else:
                        rotatablelst.append([atom1,atom2])
        except:
            pass
    #check and eliminate methyl groups
    rotatablelst = [item for item in rotatablelst if (check_ismethyl(molecule,item[0])==False and check_ismethyl(molecule,item[1])==False)]    
    for rotatable in rotatablelst:
        print "Rotatable bond: \t%i-\t%i"%(rotatable[0],rotatable[1])
    return rotatablelst 
      
def check_ismethyl(molecule, atomid):
    sumCH = 0
    if molecule.initialconformer.atomlst[atomid-1].symbol=="C":
        for i in range(0,len(molecule.connectionmatrix[atomid-1])):
            
            if molecule.connectionmatrix[atomid-1][i]=="1" and molecule.initialconformer.atomlst[i].symbol=="H":
                    sumCH+=1
    if sumCH>2:
        return True
    else:
        return False
        
    
def identify_dihedral(molecule,mask=None):
    dihedrallst = []
    matrix = molecule.connectionmatrix
    rotatablelst = identify_rotatable2(molecule,mask)
    #~ rotatablelst = identify_rotatable(molecule.filename,mask)
    # bond between atom2 -- atom3
    for rotatable in rotatablelst:
        atom2=rotatable[0]
        atom3=rotatable[1]
        #search atom1
        for i in range(0,len(matrix[atom3])):
            atom=i+1
            #if atom is not atom3 and has a connection with atom2
            if atom!=atom3 and matrix[atom2-1][i]!=None:
                atom1=atom
        #search atom4
        for i in range(0,len(matrix[atom2])):
            atom =i+1
            #if atom is not atom2 and has a connection with atom3
            if atom!=atom2 and matrix[atom3-1][i]!=None:
                atom4=atom
        #store dihedral 
        dihedrallst.append([atom1,atom2,atom3,atom4])
    return dihedrallst

def check_clashes(conformer):
    for atom1 in conformer.atomlst:
        for atom2 in conformer.atomlst:
            pos1 = [atom1.xcor,atom1.ycor, atom1.zcor]
            pos2 = [atom2.xcor,atom2.ycor, atom2.zcor]
            if (atom1!=atom2) and commonfunctions.calculate_distance(pos1,pos2)<0.4:
                return True
    return False


def conformer_to_xyz(conformer,xyznamefile):
    xyzfile=open(xyznamefile,"w")
    xyzfile.write("%i\n"%(len(conformer.atomlst)))
    if conformer.energy!=None:
        xyzfile.write("energy: %f\n"%(conformer.energy))
    else:
        xyzfile.write("\n")
        
    for atom in conformer.atomlst:
        xyzfile.write("%s %f %f %f\n"%(atom.symbol,atom.xcor,atom.ycor,atom.zcor))
    xyzfile.close()


# ===============================================================================================
# PYMOL FUNCTIONS
# ===============================================================================================

def modify_dihedral(atom1,atom2,atom3,atom4,value):
    try:
        cmd.set_dihedral(atom1,atom2,atom3,atom4,value)
    except:
        print "Error in pymol trying modified this dihedral %s %s %s %s %f"%(atom1,atom2,atom3,atom4,value)
        cmd.save("problem.pdb","all",0)
        exit()


def modify_structure(paramlst, moleculefilename):
    cmd.set('retain_order',1)
    #load structure
    pymol.cmd.load(os.path.abspath(os.curdir)+ "//" + moleculefilename)
    for param in paramlst:
        if "D" in param.paramtype:
            modify_dihedral("id " + str(param.atoms[0]),"id " + str(param.atoms[1]),"id " + str(param.atoms[2]),"id " + str(param.atoms[3]),param.value)
        #~ cmd.set('pdb_conect_all')
        xyzlst = cmd.get_model('all', 1).get_coord_list()
    pymol.cmd.delete("all")
    pymol.cmd.reinitialize
    return xyzlst
    
    
def get_dihedral(moleculefilename,atom1,atom2,atom3,atom4):
    cmd.set('retain_order',1)
    pymol.cmd.load(os.path.abspath(os.curdir)+ "//" + moleculefilename)
    value = pymol.cmd.get_dihedral("id " + str(atom1),"id " + str(atom2),"id " + str(atom3),"id " + str(atom4),0)
    pymol.cmd.delete("all")
    pymol.cmd.reinitialize
    return value    


def calculate_rmsd(conformer1,conformer2):
    conformer_to_xyz(conformer1,"tmp1.xyz")
    conformer_to_xyz(conformer2,"tmp2.xyz")
    pymol.cmd.load(os.path.abspath(os.curdir)+ "//tmp1.xyz")
    pymol.cmd.load(os.path.abspath(os.curdir)+ "//tmp2.xyz")

    rmsd = pymol.cmd.fit("tmp1","tmp2")
    pymol.cmd.delete("all")
    pymol.cmd.reinitialize
    return rmsd
