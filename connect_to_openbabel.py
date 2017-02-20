import pybel
import openbabel

def conformer_to_mmenergy(conformerlst, keywords="GAFF"):
    
    print "conformeros"
    print len(conformerlst)
    
    for conformer in conformerlst:

        xyzstr = "%i\n\n"%(len(conformer.atomlst))
        for atom in conformer.atomlst:
            xyzstr += "%s %f %f %f\n"%(atom.symbol,atom.xcor,atom.ycor,atom.zcor)
        
        
        mol = openbabel.OBMol()
        mymol = pybel.readstring("xyz", xyzstr)
        ff = openbabel.OBForceField.FindForceField(keywords)
        if ff == 0:
            print "Could not find MM forcefield using OpenBabel"
            exit()
        if ff.Setup(mymol.OBMol) == 0:
            print "Could not setup MM forcefield using OpenBabel"
            exit()
            # Calculate the energy
        conformer.energy = ff.Energy()
  
        
    return conformerlst

def is_inring(filename):
    isinringlst = []
    mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    format = conv.FormatFromExt(filename)
    conv.SetInAndOutFormats(format, format)
    conv.ReadFile(mol, filename)
    for atom in openbabel.OBMolAtomIter(mol):
        isinringlst.append(atom.IsInRing())
    return isinringlst
