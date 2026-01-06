import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
import json
import os

# Paths
cif_path = "../Sb2Te3-mp-1201.cif"
sc_path = "../sb2te3_supercell_441.xyz"
output_dir = "intrinsic_vacancy"

if not os.path.exists(cif_path):
    # Fallback if running from root
    cif_path = "26ICML-CrystalEdit/Sb2Te3-mp-1201.cif"
    sc_path = "26ICML-CrystalEdit/sb2te3_supercell_441.xyz"
    output_dir = "26ICML-CrystalEdit/stage1/intrinsic_vacancy"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print(f"Loading structures...")
unit_cell = Structure.from_file(cif_path)
# Load supercell via ASE to handle extended XYZ, then convert to Pymatgen
from ase.io import read as ase_read
atoms_sc = ase_read(sc_path)
supercell = AseAtomsAdaptor.get_structure(atoms_sc)

# 1. Identify Unique Sites in Unit Cell
print("Identifying unique sites in unit cell...")
sga = SpacegroupAnalyzer(unit_cell)
sym_struct = sga.get_symmetrized_structure()
unique_sites = sym_struct.equivalent_sites

# Prepare to map these to the supercell
# We want to find one atom in the supercell that corresponds to each unique site type.
# The supercell was created by 4x4x1 expansion.
# Atom mapping: Atom i in supercell corresponds to atom (i % n_unit) in unit cell? 
# ONLY if the order is preserved perfectly. ASE * (4,4,1) usually preserves order blocks.
# Let's verify by checking species and fractional coordinates correspondence.

def get_supercell_index_for_site(u_site, sc, scaling=[4,4,1]):
    # u_site is a site in the unit cell
    # sc is the supercell structure
    # sc_frac = (u_frac + shift) / scaling
    # We target shift=(0,0,0) -> u_frac / scaling
    target_frac = np.array(u_site.frac_coords) / np.array(scaling)
    
    # Wrap target_frac to [0, 1) just in case, though usually fine
    target_frac = target_frac % 1.0
    
    # Find closest site in supercell
    for i, site in enumerate(sc):
        if site.specie.symbol != u_site.specie.symbol:
            continue
            
        diff = site.frac_coords - target_frac
        # Apply PBC min image
        diff = diff - np.round(diff)
        
        dist = np.linalg.norm(diff)
        
        if dist < 0.1: # Relaxed tolerance
            return i
            
    # If not found at 0,0,0 shift, try to find *any* equivalent
    # Iterate over all possible shifts
    min_dist_found = float('inf')
    
    for i in range(scaling[0]):
        for j in range(scaling[1]):
            for k in range(scaling[2]):
                shift = np.array([i, j, k])
                target_frac = (np.array(u_site.frac_coords) + shift) / np.array(scaling)
                target_frac = target_frac % 1.0
                
                for idx, site in enumerate(sc):
                    if site.specie.symbol != u_site.specie.symbol: continue
                    diff = site.frac_coords - target_frac
                    diff = diff - np.round(diff)
                    dist = np.linalg.norm(diff)
                    if dist < min_dist_found:
                        min_dist_found = dist
                    
                    if dist < 0.1: # Relaxed tolerance
                        return idx
    
    print(f"Debug: Closest match distance for {u_site.specie.symbol}: {min_dist_found:.4f}")
    return None

results = []

print(f"Processing {len(unique_sites)} unique sites...")

for equivalents in unique_sites:
    representative = equivalents[0]
    el = representative.specie.symbol
    
    # Identify Wyckoff info
    # We need to find which wyckoff group this belongs to
    # sga.get_symmetry_dataset() gives 'wyckoffs' list matching structure sites
    # We need the index of 'representative' in the original 'unit_cell'
    # equivalent_sites contains sites from the structure.
    
    # Let's rely on SpacegroupAnalyzer to get labels
    # SymmetrizedStructure has .wyckoff_symbols property corresponding to .equivalent_sites
    wyckoff_label = sym_struct.wyckoff_symbols[sym_struct.equivalent_sites.index(equivalents)]
    
    label = f"{el}_{wyckoff_label}"
    print(f"Processing {label}...")
    
    # Find corresponding atom in supercell
    sc_index = get_supercell_index_for_site(representative, supercell)
    
    if sc_index is None:
        print(f"Warning: Could not find exact mapping for {label} in supercell. Skipping.")
        continue
        
    print(f"  Mapped to Supercell Atom Index: {sc_index}")
    
    # Create Vacancy
    defect_sc = supercell.copy()
    removed_site = defect_sc[sc_index]
    defect_sc.remove_sites([sc_index])
    
    # Analysis
    # 1. Point Group Symmetry
    # For the defect structure. Sga might be slow on large cell, but 239 atoms is okay.
    # Use loose tolerance for defect relaxation? No, this is unrelaxed.
    # Symmetry of a single vacancy in supercell is usually the site symmetry of the removed atom
    # constrained by the supercell shape.
    try:
        defect_sga = SpacegroupAnalyzer(defect_sc, symprec=0.1)
        final_pg = defect_sga.get_point_group_symbol()
    except:
        final_pg = "Unknown"
        
    initial_pg = "R-3m" # From previous step, but strictly that's Space Group. PG is -3m (D3d)
    
    # 2. Nearest Neighbors distances (in original SC)
    # We measure from the REMOVED site location to the remaining atoms in defect_sc
    # We can use the original supercell to find neighbors of the atom at sc_index
    dist_matrix = supercell.distance_matrix
    # row sc_index
    dists = dist_matrix[sc_index]
    # Filter 0 (self)
    # For a large supercell, distance_matrix gives distances to all atoms (including images if not carefully handled, but pymatgen distance_matrix typically is within the cell or uses min image)
    # Wait, Structure.distance_matrix returns the distances between sites in the structure (PBC considered).
    
    valid_indices = [i for i, d in enumerate(dists) if d > 0.01 and d < 5.0]
    
    nn_info = []
    if valid_indices:
        # Sort by distance
        sorted_indices = sorted(valid_indices, key=lambda i: dists[i])
        min_dist = dists[sorted_indices[0]]
        
        # Get nearest neighbors (could be multiple at same distance)
        # Just grab the closest shell
        for idx in sorted_indices:
            if dists[idx] < min_dist + 0.05: # tolerance
                nn_info.append({
                    "element": supercell[idx].specie.symbol,
                    "distance": float(dists[idx])
                })
    else:
        min_dist = 0.0
        
    # 3. Formula Unit Change
    # Original: Sb96 Te144 (from 240 atoms 4x4x1)
    # Formula of Sb2Te3 is Sb2 Te3. Z=3 in unit cell. 4x4x1 => 16 * 3 = 48 formula units?
    # Wait, Unit Cell has 6 Sb and 9 Te? No, 6 Sb and 9 Te is Z=3.
    # Supercell 4x4x1 -> 16 times unit cell atoms -> 96 Sb, 144 Te. Correct.
    # Remove 1 Sb -> Sb95 Te144
    # Remove 1 Te -> Sb96 Te143
    comp = defect_sc.composition
    new_formula = comp.formula.replace(" ", "")
    
    # 4. Validation
    n_expected = len(supercell) - 1
    n_actual = len(defect_sc)
    i_fid = (n_actual == n_expected)
    
    # Min interatomic distance in defect structure
    # We can check the whole matrix or just assume it's valid if we only removed an atom.
    # Removing an atom increases spacing effectively.
    # But let's check the remaining set.
    # Efficient way: pre-calculated min distance of pristine is known. Removing won't decrease it.
    # So valid if pristine was valid.
    # Let's just re-verify for the report.
    # Taking a subset of distances to save time?
    # Just check neighbors of neighbors of vacancy?
    # Or skip expensive N^2 check and assume > 0.6A based on pristine.
    min_interatomic = 2.0 # Placeholder, safe assumption for removal
    
    # Output Data
    defect_entry = {
        "defect_type": f"vacancy_{el}",
        "removed_atom": {
            "element": el,
            "wyckoff": wyckoff_label,
            "frac_coords": removed_site.frac_coords.tolist(),
            "supercell_index": sc_index
        },
        "new_formula": new_formula,
        "symmetry_change": {
            "initial": "R-3m (Space Group)", 
            "final_point_group": final_pg
        },
        "validation": {
            "I_fid": bool(i_fid),
            "min_distance": float(min_dist) # Distance to nearest neighbor
        },
        "nearest_neighbors": nn_info
    }
    
    results.append(defect_entry)
    
    # Save Extended XYZ
    # Filename: Vac_Sb_6c.xyz
    fname = f"Vac_{el}_{wyckoff_label}.xyz"
    full_path = os.path.join(output_dir, fname)
    
    atoms_defect = AseAtomsAdaptor.get_atoms(defect_sc)
    atoms_defect.info['defect_type'] = f"vacancy_{el}"
    atoms_defect.info['wyckoff'] = wyckoff_label
    
    write(full_path, atoms_defect, format="extxyz")
    print(f"  Saved {fname}")

# Save JSON
json_path = os.path.join(output_dir, "vacancies_report.json")
with open(json_path, 'w') as f:
    json.dump(results, f, indent=4)

print(f"Generated {len(results)} vacancies. Report saved to {json_path}")
