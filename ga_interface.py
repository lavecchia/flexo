"""ga_interface.py
connect simulation with genetic library.
"""
from objects import *
from genetic import *
import numpy as np
from commonfunctions import *
from mainfunctions import *


#~ class memberobj(object):
    #~ def __init__(self, paramlst, energy=None):
        #~ self.paramlst = paramlst
        #~ self.energy = energy
        #~ self.fitness = None
        #~ 
    #~ def set_energy(self,newenergy):
        #~ self.energy = newenergy
#~ 
    #~ def set_fitness(self,newfitness):
        #~ self.fitness = newfitness
        #~ 
        
class SetsGeneration(object):
    """Use to connect with genetic library.
    Keyword arguments:
    memberlst -- list of member objects
    minlst -- list with minima values with same format of parameters chromosome
    maxlst -- list with maxima values with same form of parameters chromosome
    """
    def __init__(self, conformerlst, molecule):
        self.conformerlst = conformerlst
        self.molecule = molecule
        minlst = []
        maxlst = []
        # take min y max values of the firts member
        for param in conformerlst[0].paramlst:
            minlst.append(param.minvalue)
            maxlst.append(param.maxvalue)

        chromosomelst = []
        fitnesslst = []
        for conformer in self.conformerlst:
            fitnesslst.append(conformer.fitness)
            chromosome = []
            for param in conformer.paramlst:
                chromosome.append(param.value)
            chromosomelst.append(chromosome)

        self.minlst = minlst 
        self.maxlst = maxlst
        self.chromosomelst = chromosomelst
        self.fitnesslst = calculate_fitness(conformerlst)
        


        
    def next(self):
        """ Return a list of dictionaries with new parameter values
        """
        ga = Generation(self.chromosomelst, self.fitnesslst,
                self.minlst, self.maxlst)
        newgeneration = ga.next()
        newmemberlst = []
        for newmember in newgeneration:
            newparamlst = []
            newvalueslst=newmember.get_chrom()
            for index in range(0,len(self.conformerlst[0].paramlst)):
                paramtype=self.conformerlst[0].paramlst[index].paramtype
                atoms=self.conformerlst[0].paramlst[index].atoms
                minvalue=self.conformerlst[0].paramlst[index].minvalue
                maxvalue=self.conformerlst[0].paramlst[index].maxvalue
                value=newvalueslst[index]
                newparamlst.append(Parameter(paramtype,atoms,value,minvalue,maxvalue))

            newconformer = self.molecule.params_to_conformer(newparamlst)
            if check_clashes(newconformer) == False:
                newmemberlst.append(newconformer)
              

                
            
        #compare input list of members (memberlst) with the new list of members (newmemberlst), to avoid calculate again the energy of new old members
        for member2 in newmemberlst:
            for member1 in self.conformerlst:
                check = 0
                for param1,param2 in map(None, member1.paramlst, member2.paramlst):
                    if param2.value == param1.value:
                        check+=1
                if check == len(member2.paramlst):
                    member2.energy = member1.energy
                    member2.atomlst = member1.atomlst
                    member1.sasa = member2.sasa 
        return newmemberlst


def calculate_fitness(memberlst):
    #fitness of J. A. Niesse and H. R. Mayne, J. Chem. Phys. https://doi.org/JCPSA6 105, 4700 (1996).
    #http://aip.scitation.org/doi/pdf/10.1063/1.472311
    # Fi = (Ei - E_min)/(E_max - E_min), fi = Fi/Ftot
    fitnesslst = []
    energylst = []
    Flst = []
    
    memberlst = [x for x in memberlst if x.energy is not None]
    
    
    for member in memberlst:
        energylst.append(member.energy)
        
    minenergy = min(energylst)
    maxenergy = max(energylst)
    
    for energy in energylst:
        Flst.append((energy-minenergy)/(maxenergy-minenergy)+0.05)
    Ftot = sum(Flst)
    
    for Fi in Flst:
        fitnesslst.append(Ftot/(Fi))
        
    for member, fitness, F in map(None, memberlst, fitnesslst, Flst):
        member.fitnesslst = fitness
        #~ print "energy: %f, fitness %f F %f"%(member.energy,fitness,F)
    return fitnesslst




    
