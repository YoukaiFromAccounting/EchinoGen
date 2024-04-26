from pymol import cmd

def load_and_label(file_path, residue_numbers):
    strut_name = load_in(file_path)
    label_all(strut_name, residue_numbers)

# load in the pdb file with file path
def load_in(file_path):
    strut_name = file_path.split('/')[-1].split('.')[0]
    cmd.load(file_path, strut_name)
    return strut_name

# Select and label residues
def label_all(object_name, residue_numbers):
    residue_numbers = residue_numbers.split(' ')
    residue_numbers = '+'.join(residue_numbers)
    selection_command = f"resi {residue_numbers} and object {object_name}"
    cmd.select("residues", selection_command)
    cmd.color("red", "residues")
    cmd.show("licorice", "residues")


# Extend cmd with new function
cmd.extend("load_and_label", load_and_label)

# Setup command-line argument auto-completion
cmd.auto_arg[1]['load_and_label'] = [lambda: cmd.Shortcut("../data/idoform1.pdb"), 'file_path', '']
cmd.auto_arg[2]['load_and_label'] = [lambda: cmd.Shortcut("12 13 14"), 'residue_numbers', '']