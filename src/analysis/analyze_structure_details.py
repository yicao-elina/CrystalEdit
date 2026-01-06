import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.io.cif import CifWriter
import os

# Paths
cif_path = "data/structures/Sb2Te3-mp-1201.cif"
xyz_path = "data/structures/sb2te3_supercell_441.xyz"

if not os.path.exists(cif_path):
    cif_path = "Sb2Te3-mp-1201.cif"
if not os.path.exists(xyz_path):
    xyz_path = "sb2te3_supercell_441.xyz"

print(f"--- Structure Analysis Report ---")

# 1. Load Unit Cell (for Symmetry/Wyckoff)
print(f"\nLoading Unit Cell from {cif_path}...")
structure = Structure.from_file(cif_path)

# 2. Symmetry Analysis
sga = SpacegroupAnalyzer(structure)
dataset = sga.get_symmetry_dataset()

sg_symbol = dataset["international"]
sg_number = dataset["number"]
print(f"\nSpace Group: {sg_symbol} ({sg_number})")

lat = structure.lattice
print(f"\nLattice Parameters:")
print(f"a = {lat.a:.4f} Å")
print(f"b = {lat.b:.4f} Å")
print(f"c = {lat.c:.4f} Å")
print(f"α = {lat.alpha:.2f}°")
print(f"β = {lat.beta:.2f}°")
print(f"γ = {lat.gamma:.2f}°")

# 3. Wyckoff Positions & Unique Sites
print(f"\nSymmetrically Unique Atomic Sites & Wyckoff Positions:")
print(f"{ 'Element':<10} { 'Wyckoff':<10} { 'Multiplicity':<15} { 'Fract Coords (x, y, z)':<30}")
print("-" * 70)

# Pymatgen equivalent_sites returns list of list of sites
# We need to map them to Wyckoff labels.
# sga.get_symmetrized_structure() gives a SymmetrizedStructure object
sym_struct = sga.get_symmetrized_structure()

# For CrystalNN
cnn = CrystalNN()

unique_sites_info = []

for i, equivalent_sites in enumerate(sym_struct.equivalent_sites):
    site = equivalent_sites[0] # Representative
    el = site.specie.symbol
    wyckoff = sym_struct.wyckoff_symbols[i]
    multiplicity = len(equivalent_sites)
    coords = f"[{site.frac_coords[0]:.4f}, {site.frac_coords[1]:.4f}, {site.frac_coords[2]:.4f}]"
    
    print(f"{el:<10} {wyckoff:<10} {multiplicity:<15} {coords:<30}")
    
    # Coordination Analysis
    # Analyze representative site
    cn_info = cnn.get_nn_info(structure, structure.index(site))
    cn = len(cn_info)
    neighbors = [n['site'].specie.symbol for n in cn_info]
    neighbor_summary = ", ".join([f"{n}{neighbors.count(n)}" for n in sorted(set(neighbors))])
    
    unique_sites_info.append({
        "Element": el,
        "Wyckoff": wyckoff,
        "CN": cn,
        "Polyhedron": neighbor_summary
    })

# 4. Coordination Polyhedra Analysis
print(f"\nCoordination Polyhedra Analysis (CN):")
for info in unique_sites_info:
    print(f"  {info['Element']} (Wyckoff {info['Wyckoff']}): CN = {info['CN']} (Neighbors: {info['Polyhedron']})")

# 5. Original CIF Content
print(f"\nOriginal CIF Metadata (First 20 lines):")
with open(cif_path, 'r') as f:
    lines = f.readlines()
    for line in lines[:20]:
        print(line.strip())

# 6. Verify Supercell
print(f"\nSupercell Analysis ({xyz_path}):")
try:
    # Use ASE to read extended XYZ lattice info
    from ase.io import read
    atoms = read(xyz_path)
    print(f"  Total Atoms: {len(atoms)}")
    print(f"  Chemical Formula: {atoms.get_chemical_formula()}")
    cell = atoms.get_cell()
    print(f"  Supercell Lattice: {cell.cellpar()}")
except Exception as e:
    print(f"  Could not read supercell: {e}")
