import numpy as np
import os

# Handle numpy 2.0 compatibility for older ASE versions
if not hasattr(np, 'product'):
    np.product = np.prod

from ase.io import read, write

# Set file paths
cif_file = "Sb2Te3-mp-1201.cif"
output_file = "sb2te3_supercell_441.xyz"

# 1. Read the pristine structure
print(f"Reading {cif_file}...")
atoms = read(cif_file)

# 2. Build 4x4x1 supercell
# Multiplying atoms object by (4, 4, 1) repeats the unit cell
print("Building 4x4x1 supercell...")
supercell = atoms * (4, 4, 1)

# 3. Save in extended xyz format
print(f"Saving supercell to {output_file}...")
write(output_file, supercell, format="extxyz")

print("Done.")
