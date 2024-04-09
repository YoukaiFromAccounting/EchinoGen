#load library
from pymol import cmd, stored
import sys
import os

#function to judge if a file exists
def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

#function to switch atom names within a residue
def switchName(residueSelection, atomName1, atomName2):
  """
  switch the name of two atoms
  """
  cmd.alter(residueSelection+" and name "+atomName1, 'name="gzt"')
  cmd.alter(residueSelection+" and name "+atomName2, 'name="'+atomName1+'"')
  cmd.alter(residueSelection+" and name gzt", 'name="'+atomName2+'"')


#function to change atom names of some residues with symetric sidechain
def flipAtomName(targetResidueSelection):
  """
  switch the atom names of specific residues
  """
# Create flipped residue
  cmd.create("flippedRes",targetResidueSelection+" and not alt B")
  targetResidueCa=cmd.get_model("flippedRes and name CA")
  for g in targetResidueCa.atom:
#   print g.resn
    if g.resn=='ARG':
     switchName("flippedRes", "NH1", "NH2")
    elif g.resn=='HIS':
     switchName("flippedRes", "ND1", "CD2")
     switchName("flippedRes", "CE1", "NE2")
    elif g.resn=='ASP':
     switchName("flippedRes", "OD1", "OD2")
    elif g.resn=='PHE':
     switchName("flippedRes", "CD1", "CD2")
     switchName("flippedRes", "CE1", "CE2")
    elif g.resn=='GLN':
     switchName("flippedRes", "OE1", "NE2")
    elif g.resn=='GLU':
     switchName("flippedRes", "OE1", "OE2")
    elif g.resn=='LEU':
     switchName("flippedRes", "CD1", "CD2")
    elif g.resn=='ASN':
     switchName("flippedRes", "OD1", "ND2")
    elif g.resn=='TYR':
     switchName("flippedRes", "CD1", "CD2")
     switchName("flippedRes", "CE1", "CE2")
    elif g.resn=='VAL':
      switchName("flippedRes", "CG1", "CG2")
  cmd.sort()
# cmd.label("flippedRes","name")
  return "flippedRes"

#main function
def rmsdByRes(referenceProteinChain, sel, targetProteinChain):
  """
  Update
    Zhenting Gao on 7/28/2016

  USAGE

    rmsf referenceProteinChain, targetProteinChain, selection [,byres=0], [reference_state=1]

    Calculate the RMSD for each residue pairs from two chains of the same protein from two crystal structures.

  Workflow
    Read reference and target pdb files
    Align two structures
        sel target, proA and chain A
            #define target protein chain
        sel refrence, proB and chain A
            #define reference protein chain
        align target, reference
            #automatical alignment
    Clean attributes
        otherwise rms_cur will fail


  """
# Create temporary objects, exclude alternative conformation B
  cmd.create("ref_gzt", referenceProteinChain+" and polymer and not alt B")
  cmd.alter("ref_gzt", "chain='A'")
  cmd.alter("ref_gzt", "segi=''")
  cmd.create("target_gzt", targetProteinChain+" and polymer and not alt B")
  cmd.alter("target_gzt", "chain='A'")
  cmd.alter("target_gzt", "segi=''")
#  cmd.align("target_gzt","ref_gzt",object="align")
# parameters
  outputText=""
  res2Check=['HIS','ASP','ARG','PHE','GLN','GLU','LEU','ASN','TYR','VAL']
       
# select alpha carbon of selected residues in reference structure
  calpha=cmd.get_model(sel+" and name CA and not alt B")
  
  for g in calpha.atom:
#  print g.resi+g.resn
   if cmd.count_atoms("ref_gzt and polymer and resi "+g.resi)==cmd.count_atoms("target_gzt and polymer and resi "+g.resi):
    rmsdRes=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi,"target_gzt and polymer and resi "+g.resi)
    rmsdResCa=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi+" and name ca","target_gzt and polymer and resi "+g.resi+" and name ca")
    rmsdResBackbone=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi+" and name ca+n+c+o","target_gzt and polymer and resi "+g.resi+" and name ca+n+c+o")
#  calculate minimum rmsd
    rmsdResMin=rmsdRes
    if g.resn in res2Check:
      flippedRes=flipAtomName("target_gzt and polymer and resi "+g.resi)
      rmsdFlippedRes=cmd.rms_cur("ref_gzt and polymer and resi "+g.resi,flippedRes)
      if rmsdFlippedRes<rmsdRes:
       rmsdResMin=rmsdFlippedRes

#    print cmd.count_atoms("ref_gzt and polymer and resi "+g.resi),cmd.count_atoms("target_gzt and polymer and resi "+g.resi)
    outputText+="%s,%s,%s,%.3f,%.3f,%.3f,%.3f\n" % (targetProteinChain,g.resn,g.resi,rmsdRes,rmsdResCa,rmsdResBackbone,rmsdResMin)
  
  print(outputText)

# Destroy temporary objects
  cmd.delete("ref_gzt target_gzt align res_gzt "+flippedRes)
  
# Save data into csv
  outputFile='rmsdByRes_'+sel+'.csv'
  f=open(outputFile,'a')
  if not is_non_zero_file(outputFile):
   f.write("targe,residueName,residueId,allAtomRMSD,rmsdResCa,rmsdResBackbone,allAtomRMSDMin\n")
  f.write(outputText)
  f.close()
  print("Results saved in "+outputFile)
  
cmd.extend("rmsdByRes",rmsdByRes)