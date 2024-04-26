from pymol import cmd

# Select and label residues
def label_all(object_name, residue_numbers):
    residue_numbers = residue_numbers.split(' ')
    residue_numbers = '+'.join(residue_numbers)
    selection_command = f"resi {residue_numbers} and object {object_name}"
    cmd.select("residues", selection_command)
    cmd.color("red", "residues")
    cmd.show("licorice", "residues")


# Extend cmd with new function
cmd.extend("label_all", label_all)

# Setup command-line argument auto-completion
cmd.auto_arg[1]['label_all'] = [cmd.object_sc, 'object_name', '']
cmd.auto_arg[2]['label_all'] = [lambda: cmd.Shortcut("12 13 14"), 'residue_numbers', '']