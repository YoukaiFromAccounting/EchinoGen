import sys
import pymol
from pymol import cmd

# Initialize PyMOL if not already running in a script
if not hasattr(sys, 'argv'):
    pymol.finish_launching()

# Function to load a file and label residues based on command arguments
def load_and_label(file_name, residue_numbers):
    # Parse residue numbers assuming they are comma-separated
    residue_numbers = residue_numbers.split(' ')
    structure_name = load_in_pdb(file_name)
    label_residues(structure_name, residue_numbers)

# Load the PDB file
def load_in_pdb(file_name):
    parts = file_name.split('.')
    structure_name = parts[0] + "_structure"
    cmd.load(file_name, structure_name)
    return structure_name

# Select and label residues
def label_residues(structure_name, residue_numbers):
    residue_numbers_str = '+'.join(residue_numbers)
    cmd.select("residues", f"{structure_name} and resi {residue_numbers_str}")
    cmd.label("residues", "name+resi")

# Extend cmd with new function
cmd.extend("load_and_label", load_and_label)
cmd.load_and_label = load_and_label

# Setup command-line argument auto-completion
cmd.auto_arg[1]['load_and_label'] = [lambda: cmd.Shortcut("isoform1.pdb"), 'file_name', '']
cmd.auto_arg[2]['load_and_label'] = [lambda: cmd.Shortcut("12 13 14"), 'residue_numbers', '']