from prody import *
from Bio.PDB import PDBParser, Superimposer
import numpy as np
import matplotlib.pyplot as plt
import sys

# Accept file names as args
file1 = sys.argv[1]
file2 = sys.argv[2]

# Parse PDB files and align sequences
iso1 = parsePDB(file1)
iso2 = parsePDB(file2)
matches = matchChains(iso1, iso2)

# Isolate matched sequenced 
iso1_match = matches[0][0] 
iso2_match = matches[0][1] 

# Write new pdb with only the aligned portions for comparison
writePDB('new_isoform1.pdb',iso1_match)
writePDB('new_isoform2.pdb',iso2_match)

def calculate_rmsd(structure1, structure2):
    si = Superimposer()
    #extract CA coordinates
    fixed_atoms = extract_coordinates(structure1)
    moving_atoms = extract_coordinates(structure2)

    #check that the two structures have the same number of atoms
    if len(fixed_atoms) != len(moving_atoms):
        raise ValueError("Structures do not have the same number of CA atoms")

    #calculate the RMSD
    fixed_coords = np.array([atom.coord for atom in fixed_atoms])
    moving_coords = np.array([atom.coord for atom in moving_atoms])
    si.set_atoms(fixed_atoms, moving_atoms)
    si.apply(moving_atoms)
    total_rmsd = si.rms

    #this command will calculate the residue wise RMSD
    residue_rmsds = np.sqrt(np.sum((fixed_coords - moving_coords)**2, axis=1))
    return total_rmsd, residue_rmsds

def extract_coordinates(structure):
    #extract the C alpha coordinates for backbone comparison
    atoms = []
    for residue in structure.get_residues():
        if 'CA' in residue:
            atoms.append(residue['CA'])
    return atoms

def read_structure(filename):
    parser = PDBParser()
    structure = parser.get_structure(id(filename), filename)
    return structure

iso_1_structure = read_structure('new_isoform1.pdb')
iso_2_structure = read_structure('new_isoform2.pdb')

# Calculate RMSDs
total_rmsd, residue_rmsds = calculate_rmsd(iso_1_structure, iso_2_structure)

indices = []
RMSD_vals = []
for i, rmsd in enumerate(residue_rmsds, 1):
    indices.append(i)
    RMSD_vals.append(rmsd)

# Plotting
plt.plot(indices, RMSD_vals)

# Extract file names for plot label 
file1label1 = str(file1)
file1label_final = file1label1[:-4]
file2label1 = str(file2)
file2label_final = file2label1[:-4]

# Adding labels and title
plt.xlabel('Index')
plt.ylabel('RMSD')
plt.title(file1label_final+' vs '+file2label_final +' RMSD Comparison')
plt.legend()  # Show legend

# Save plot as PNG
plt.savefig('RMSD_residue_plot_'+file1label1+'_'+file2label1+'.png')