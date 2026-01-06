import numpy as np
import json
import os
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write as ase_write
from pymatgen.analysis.local_env import CrystalNN

# Paths
cif_path = "data/structures/Sb2Te3-mp-1201.cif"
output_dir = "data/stage2/surface_slabs"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print("Loading unit cell...")
unit_cell = Structure.from_file(cif_path)

# Add oxidation states for analysis
unit_cell.add_oxidation_state_by_element({"Sb": 3, "Te": -2})

# 1. Slab Generation
print("Generating (001) slabs...")
# We want 5 QLs.
# The unit cell (conventional) has 3 QLs and c ~ 30.5 A.
# 1 QL ~ 10.1 A.
# 5 QLs ~ 50.5 A.
# We set min_slab_size to cover 5 QLs.
min_slab_size = 48.0 # Angstroms

slab_gen = SlabGenerator(
    unit_cell, 
    miller_index=(0, 0, 1), 
    min_slab_size=min_slab_size, 
    min_vacuum_size=15.0, 
    center_slab=True,
    reorient_lattice=True 
)

slabs = slab_gen.get_slabs()
print(f"Generated {len(slabs)} unique slab terminations.")

# Identify Te-terminated and Sb-terminated
te_slab = None
sb_slab = None

def analyze_termination(slab):
    cart_coords = slab.cart_coords
    z = cart_coords[:, 2]
    max_z_idx = np.argmax(z)
    min_z_idx = np.argmin(z)
    
    top_el = slab[max_z_idx].specie.symbol
    bottom_el = slab[min_z_idx].specie.symbol
    
    return top_el, bottom_el

def count_qls(slab):
    num_atoms = len(slab)
    # QL count = num_atoms / 5 (if stoichiometric 1x1).
    return num_atoms / 5.0

final_slabs = []

for i, slab in enumerate(slabs):
    top, bottom = analyze_termination(slab)
    ql_count = count_qls(slab)
    print(f"Slab {i}: {top}-{bottom}, {ql_count:.1f} QLs, Formula: {slab.composition.formula}")
    
    # 1. Te-terminated (Symmetric)
    # Expected: Te top, Te bottom.
    if top == "Te" and bottom == "Te":
        # Check stoichiometry (Sb2Te3 reduced is Sb2Te3)
        if slab.composition.reduced_formula == unit_cell.composition.reduced_formula:
             # Unit cell x2 was 6 QL. Slab might be 6 QL.
             if abs(ql_count - 6.0) < 0.5 or abs(ql_count - 5.0) < 0.5:
                 te_slab = slab
                 print(f"  -> Identified Te-terminated slab.")

    # 2. Sb-terminated
    if top == "Sb" and bottom == "Sb":
        sb_slab = slab
        print(f"  -> Identified Sb-terminated slab.")

# Fallback: Pick Slab 2 (Te-Te) if not identified
if te_slab is None:
    candidates = [s for s in slabs if analyze_termination(s) == ("Te", "Te")]
    if candidates:
        te_slab = candidates[0]
        print("Picked fallback Te-Te slab.")

# Ensure 5 QL thickness for Te slab
if te_slab is not None:
    # Scale to supercell first
    scaling = [4, 4, 1]
    te_slab.make_supercell(scaling)
    print("Scaled Te-slab to 4x4.")
    
    # Check thickness
    # 1 QL = 5 atoms * 16 (4x4) = 80 atoms.
    # Current slab might be 6 QL (480 atoms).
    # Target: 5 QL (400 atoms).
    n_atoms = len(te_slab)
    target_atoms = 5 * 5 * 16 # 5 atoms/QL * 16 per layer? No.
    # Unit cell has 5 atoms per QL (Sb2Te3).
    # Supercell 4x4 has 16 formula units per layer?
    # 4x4x1 supercell of (Sb2Te3) => 16 * 5 = 80 atoms per QL.
    # 5 QLs = 400 atoms.
    
    if n_atoms > 400:
        to_remove_count = n_atoms - 400
        if to_remove_count > 0 and to_remove_count % 80 == 0:
            print(f"Trimming {to_remove_count} atoms to get 5 QL...")
            # Sort by z
            sorted_sites = sorted(enumerate(te_slab), key=lambda k: k[1].coords[2])
            # Remove from top (last indices in sorted list)
            to_remove_indices = [x[0] for x in sorted_sites[-to_remove_count:]]
            te_slab.remove_sites(to_remove_indices)
            print("Trimmed to 5 QL.")
        else:
            print(f"Warning: Slab size {n_atoms} not multiple of 80 or not easily trimmed to 400.")

if sb_slab is None and te_slab is not None:
    print("Constructing Sb-terminated slab from Te-terminated...")
    sb_slab = te_slab.copy()
    
    # Identify highest and lowest z atoms
    cart_coords = sb_slab.cart_coords
    z = cart_coords[:, 2]
    z_max = np.max(z)
    z_min = np.min(z)
    
    # Remove top Te layer and bottom Te layer
    # Te1 layer is the boundary.
    # In 4x4, 16 atoms per layer.
    to_remove = []
    for idx, site in enumerate(sb_slab):
        if site.coords[2] > z_max - 1.0 or site.coords[2] < z_min + 1.0:
            if site.specie.symbol == "Te":
                to_remove.append(idx)
    
    sb_slab.remove_sites(to_remove)
    print(f"Removed {len(to_remove)} Te atoms to expose Sb.")
    
elif sb_slab is not None:
    sb_slab.make_supercell([4, 4, 1])

# Save Files
if te_slab:
    ase_te = AseAtomsAdaptor.get_atoms(te_slab)
    ase_te.info['termination'] = "Te"
    ase_te.info['layers'] = "5 QL"
    ase_write(os.path.join(output_dir, "Sb2Te3_001_Te_Terminated_5QL.xyz"), ase_te, format="extxyz")
    te_slab.to(filename=os.path.join(output_dir, "Sb2Te3_001_Te_Terminated_5QL.cif"))

if sb_slab:
    ase_sb = AseAtomsAdaptor.get_atoms(sb_slab)
    ase_sb.info['termination'] = "Sb"
    ase_sb.info['layers'] = "Sb-terminated"
    ase_write(os.path.join(output_dir, "Sb2Te3_001_Sb_Terminated.xyz"), ase_sb, format="extxyz")
    sb_slab.to(filename=os.path.join(output_dir, "Sb2Te3_001_Sb_Terminated.cif"))

print("Slabs generated and saved.")

# 2. Analysis
print("Performing Surface Analysis...")

results = {
    "Te_terminated": {},
    "Sb_terminated": {}
}

def analyze_slab_properties(slab, label):
    if slab is None: return
    
    # Sort by z
    sorted_sites = sorted(enumerate(slab), key=lambda x: x[1].coords[2])
    z_coords = [s.coords[2] for s in slab]
    z_max = max(z_coords)
    z_min = min(z_coords)
    
    surface_atoms_data = []
    cnn = CrystalNN()
    
    for i, site in sorted_sites:
        z = site.coords[2]
        is_surf = False
        if z > z_max - 2.5:
            is_surf = True
            loc = "Top"
        elif z < z_min + 2.5:
            is_surf = True
            loc = "Bottom"
            
        if is_surf:
            try:
                nn_info = cnn.get_nn_info(slab, i)
                cn = len(nn_info)
            except:
                cn = 0
            
            el = site.specie.symbol
            expected_cn = 6 if el == "Sb" else 3
            deficit = max(0, expected_cn - cn)
            
            surface_atoms_data.append({
                "index": i,
                "element": el,
                "location": loc,
                "cn": cn,
                "expected_cn_bulk": expected_cn,
                "cn_deficit": deficit
            })
            
    matrix = slab.lattice.matrix
    area = np.linalg.norm(np.cross(matrix[0], matrix[1]))
    
    total_deficit = sum([x['cn_deficit'] for x in surface_atoms_data])
    broken_bonds_density = (total_deficit / 2.0) / area
    
    results[label] = {
        "surface_area_A2": area,
        "broken_bond_density_per_A2": broken_bonds_density,
        "avg_cn_deficit_surface": np.mean([x['cn_deficit'] for x in surface_atoms_data]) if surface_atoms_data else 0,
        "surface_atoms_sample": surface_atoms_data[:5]
    }
    
    # Layer Spacings
    z_unique = np.sort(np.unique(np.round(z_coords, 1)))
    spacings = np.diff(z_unique)
    results[label]["layer_spacings"] = spacings.tolist()
    
    return results[label]

print("\nAnalyzing Te-terminated slab...")
analyze_slab_properties(te_slab, "Te_terminated")

print("\nAnalyzing Sb-terminated slab...")
analyze_slab_properties(sb_slab, "Sb_terminated")

with open(os.path.join(output_dir, "surface_analysis.json"), "w") as f:
    json.dump(results, f, indent=4)

print("Analysis complete.")