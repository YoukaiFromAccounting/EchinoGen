import torch
import esm
import sys
import os

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 2:
    print("Usage: python3 script.py <sequence_id>")
    sys.exit(1)

# Get the input FASTA file name from the command-line arguments
fasta_file = sys.argv[1]

# Check if the FASTA file exists
if not os.path.isfile(fasta_file):
    print("Error: File not found:", fasta_file)
    sys.exit(1)

# Load the protein sequence from the FASTA file
with open(fasta_file, "r") as f:
    # Skip the header line starting with '>'
    next(f)
    sequence = "".join(line.strip() for line in f)
    
# Generate the output PDB file name using the input FASTA file name
output_file = os.path.splitext(os.path.basename(fasta_file))[0] + ".pdb"


print("Building model...")
model = esm.pretrained.esmfold_v1()
print("Evaluating model...")
model = model.eval().cuda()

# Optionally, uncomment to set a chunk size for axial attention. This can help reduce memory.
# Lower sizes will have lower memory requirements at the cost of increased speed.
# model.set_chunk_size(128)

print("Model ready. Inferring sequence...")
with torch.no_grad():
    output = model.infer_pdb(sequence)

with open(output_file, "w") as f:
    f.write(output)

import biotite.structure.io as bsio
struct = bsio.load_structure(output_file, extra_fields=["b_factor"])
print(struct.b_factor.mean())  # this will be the pLDDT
# 88.3



"""
#Run visualization using ESM
	interact --gpu
	conda activate esmfold
	mkdir esm_results
	cd esm_results
	#Obtain the results images
	python3 ../esm_qstart.py
	#Obtain the visualization PDB file
	python3 ../esm_nstart.py
	cd ..
	
    """
    
