import pymol
from pymol import cmd, stored

# initialize PyMOL
pymol.finish_launching()

# load the PDB file
cmd.load("./isoform3.pdb", "isoform3_structure") # file name, object name

# # run the script 'findseq.py' with arguments
# cmd.run("./findseq.py")
# cmd.do("findseq EDTYY, 2wtt")

# # get the names of all loaded molecular objects
# object_names = cmd.get_names(type='selections')
# print(object_names)
# # color the selected atoms
# cmd.color("red", object_names[0])
# # change the representation of the selected atoms
# cmd.show("licorice", object_names[0])

# label the selected atoms
label_index = [33, 46, 50, 83]
for i in label_index:
    cmd.color("red", f"resi {i}")
    cmd.show("locorice", f"resi {i}")