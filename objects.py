from commonparameters import *
import commonfunctions
import copy
import mainfunctions
import connect_to_openbabel

class myAtom(object):
    def __init__(self,symbol,number,x,y,z,atomtype=None,charge=None,asa=None,radius=None,freeze=False):
        self.symbol = symbol
        self.atomnumber = number
        self.xcor = x
        self.ycor = y
        self.zcor = z
        self.asa = asa
        self.atomtype = atomtype
        self.charge = charge
        self.radius = radius
        self.freeze = freeze


class myConformer(object):
    def __init__(self,atomlst=[],paramlst=[],energy=None,sasa=None):
        self.atomlst = atomlst
        self.sasa = sasa
        self.energy = energy
        self.paramlst = paramlst    
        self.fitness = None  
    
    def set_energy(self,energy):
        self.energy=float(energy)

        
                
class myMolecule(object):
    def __init__(self,filename=None,paramlst=None):
        self.filename = filename
        self.initialparamlst = paramlst
        self.initialconformer = myConformer()
        self.bondlst = []
        self.bestconformer = bestConformerList()
        self.conformerlst = []
        self.rotatablelst = [] 
        self.isinringlst = [] 
        self.read_mol2(filename)
        self.connectionmatrix=self.calculate_connectionmatrix()
       

    def read_mol2(self,filename):
        moleculefile = open(filename,"r")
        tag = "intro"
        self.initialconformer.atomlst = []
        self.bondlst = []
        for line in moleculefile:
            if "@<TRIPOS>ATOM" in line:
                tag = "coordinates"
            elif "@<TRIPOS>BOND" in line:
                tag = "bonds"
            elif tag=="coordinates":
                linediv=line.split()
                atomnumber = int(linediv[0])
                symbol = linediv[5].split(".")[0]
                x = float(linediv[2])
                y = float(linediv[3])
                z = float(linediv[4])
                atomtype = linediv[5]
                charge= float(linediv[8])
                self.initialconformer.atomlst.append(myAtom(symbol,atomnumber,x,y,z,atomtype,charge))
            elif tag=="bonds":
                try:
                    linediv=line.split()
                    bondnumber=int(linediv[0])
                    atom1=int(linediv[1])
                    atom2=int(linediv[2])
                    bondorder=linediv[3]
                    self.bondlst.append([atom1,atom2,bondorder])
                except:
                    tag="other"
        moleculefile.close()
        self.isinringlst = connect_to_openbabel.is_inring(filename)
        return 0
            
    def params_to_conformer(self,paramlst):
        xyzlst = mainfunctions.modify_structure(paramlst,self.filename)
        atomlst=[]
        for atom,xyz in map(None,self.initialconformer.atomlst,xyzlst):
            atomlst.append(copy.copy(atom)) #problem memory
            #~ atomlst.append(atom)
            atomlst[-1].xcor=xyz[0]
            atomlst[-1].ycor=xyz[1]
            atomlst[-1].zcor=xyz[2]
        return myConformer(atomlst,paramlst)
    
    def calculate_connectionmatrix(self):
        matrixsize = len(self.initialconformer.atomlst)
        matrix = [None] * matrixsize
        for i in range(matrixsize):
            matrix[i] = [None] *  matrixsize
        # connection [0]:atom1, [1]:atom2, [2]:bond order
        for connection in self.bondlst:
            # fill matrix elements atom M,atom N and atomN,atomM with bond order
            matrix[connection[0]-1][connection[1]-1] = connection[2]
            matrix[connection[1]-1][connection[0]-1] = connection[2]
        return matrix

                        
class bestConformerList(object):
    def __init__(self,bestconfnum=BESTCONFNUM,energythreshold=ENERGYTHRESHOLD,rmsdthreshold=RMSDTHRESHOLD):
        self.bestconfnum = bestconfnum
        self.bestconformerlst = []
        self.energythreshold=energythreshold
        self.rmsdthreshold=rmsdthreshold
    
    # update the list of the best conformations
    def update_list(self,newconformerlst):
        joinlst = self.bestconformerlst + newconformerlst
        sortedconformerlst = sorted(joinlst, key=lambda x: x.energy)
        selectedconformerlst = []
        conformernumber = len(sortedconformerlst)
        for i in range(0,conformernumber):
            nselected = len(selectedconformerlst)
            #if the energy is empty, skip this new conformer
            if sortedconformerlst[i].energy == None:
                pass
            #if there is not conformer selected
            elif nselected == 0:
                selectedconformerlst.append(sortedconformerlst[i])
            elif type(self.bestconfnum) == int:
                if (0 < nselected <=self.bestconfnum):
                    #if energy difference from the last selected conformation is below energythreshold
                    if (selectedconformerlst[-1].energy <= sortedconformerlst[i].energy <= selectedconformerlst[-1].energy + self.energythreshold):
                        #then if the rmsd is lager than threshold the conformer is selected
                        acceptable = True
                        for acceptedconformer in selectedconformerlst:
                            if acceptable==False:
                                pass
                            elif mainfunctions.calculate_rmsd(acceptedconformer, sortedconformerlst[i])< self.rmsdthreshold:
                                acceptable = False
                        if acceptable==True:
                            selectedconformerlst.append(sortedconformerlst[i])
                            
                    #or if the energy difference from the last selected conformation is larger than energythreshold
                    elif (sortedconformerlst[i].energy > (selectedconformerlst[-1].energy + self.energythreshold)): 
                        selectedconformerlst.append(sortedconformerlst[i])
            
            #take conformer in a range of 2.5kcal mol-1 from most stable
            elif self.bestconfnum == "automatic" and nselected > 0 and (sortedconformerlst[i].energy <(selectedconformerlst[0].energy + 2.5)):
                if (selectedconformerlst[-1].energy <= sortedconformerlst[i].energy <= selectedconformerlst[-1].energy + self.energythreshold):
                    #then if the rmsd is lager than threshold the conformer is selected
                    rmsd = mainfunctions.calculate_rmsd(selectedconformerlst[-1], sortedconformerlst[i])
                    if rmsd >self.rmsdthreshold: 
                        print "rmsd %f"%(rmsd)
                        selectedconformerlst.append(sortedconformerlst[i])
                elif (sortedconformerlst[i].energy > (selectedconformerlst[-1].energy + self.energythreshold)): 
                    selectedconformerlst.append(sortedconformerlst[i])
        self.bestconformerlst = selectedconformerlst
        return selectedconformerlst
           
            
class Parameter(object):
    def __init__(self,paramtype,atoms,value,minvalue,maxvalue):
        self.paramtype = paramtype
        self.atoms = atoms
        self.value = value
        self.minvalue = minvalue
        self.maxvalue = maxvalue
        
    
   

