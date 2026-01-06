import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read as ase_read
from ase.io import write as ase_write
import json
import os
import random

# Paths
cif_path = "../Sb2Te3-mp-1201.cif"
sc_path = "../sb2te3_supercell_441.xyz"
output_dir = "substitutional_cr_on_sb"

if not os.path.exists(cif_path):
    # Fallback
    cif_path = "26ICML-CrystalEdit/Sb2Te3-mp-1201.cif"
    sc_path = "26ICML-CrystalEdit/sb2te3_supercell_441.xyz"
    output_dir = "26ICML-CrystalEdit/stage1/substitutional_cr_on_sb"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print("Loading structures...")
unit_cell = Structure.from_file(cif_path)
atoms_sc = ase_read(sc_path)
supercell = AseAtomsAdaptor.get_structure(atoms_sc)

# Add oxidation states for Ewald sum (Sb+3, Te-2, Cr+3)
supercell.add_oxidation_state_by_element({"Sb": 3, "Te": -2})

# 1. Symmetry Analysis of Supercell
print("Analyzing symmetry orbits in supercell...")
sga_sc = SpacegroupAnalyzer(supercell, symprec=0.01)
sym_struct_sc = sga_sc.get_symmetrized_structure()
equivalent_indices = sym_struct_sc.equivalent_indices

# Identify Sb orbits
sb_orbits = []
for orbit in equivalent_indices:
    rep_idx = orbit[0]
    site = supercell[rep_idx]
    if site.specie.symbol == "Sb":
        # Identify Wyckoff label from unit cell mapping if possible, or just use supercell index
        # In primitive R-3m Sb2Te3, Sb is at 6c.
        # In supercell, all Sb should ideally be equivalent if symmetry is preserved.
        # Let's check multiplicity.
        
        orbit_label = f"Sb_Orbit_{rep_idx}"
        
        # Determine Wyckoff
        # In unit cell Sb is 6c.
        wyckoff_label = "6c" # Default for Sb in Sb2Te3
        
        sb_orbits.append({
            "label": orbit_label,
            "wyckoff": wyckoff_label,
            "indices": orbit
        })

print(f"Identified {len(sb_orbits)} unique Sb orbits in supercell.")

# 2. Combinatorial Generation
defect_concentration = 0.05
n_atoms = len(supercell)
# Concentration calculation:
# Target is x=0.05 in Sb_{2-x}Cr_{x}Te_3 ? Or 5% of cation sites?
# Formula unit Sb2Te3 -> 5 atoms.
# Usually x is per formula unit.
# Sb_{1.95}Cr_{0.05}Te_3
# Total cation sites = 2. Total sites = 5.
# x=0.05 means 0.05/2 = 2.5% of Sb sites are replaced.
# In supercell: 96 Sb atoms.
# 2.5% of 96 = 2.4 atoms. So ~2 or 3 atoms.
# Let's generate:
# 1. Single substitution (Dilute limit)
# 2. Multiple substitution (Target x=0.05 => ~2-3 atoms)

generated_configs = []

# A. Single Substitution
print("Generating single Sb->Cr substitutions...")
for orbit in sb_orbits:
    idx = orbit['indices'][0]
    label = orbit['label']
    wyckoff = orbit['wyckoff']
    
    defect_sc = supercell.copy()
    # Substitute
    defect_sc.replace(idx, "Cr")
    
    # Analyze Local Environment
    # 1. Cr-Te distances
    # Distance matrix
    dists = defect_sc.distance_matrix[idx]
    # Filter for Te neighbors (indices of Te)
    te_indices = [i for i, s in enumerate(defect_sc) if s.specie.symbol == "Te"]
    
    cr_te_dists = []
    for te_idx in te_indices:
        d = dists[te_idx]
        if d < 3.5: # First shell roughly
            cr_te_dists.append(d)
    
    cr_te_dists.sort()
    min_dist = cr_te_dists[0] if cr_te_dists else 0.0
    nn_summary = f"CN={len(cr_te_dists)} range=[{min_dist:.3f}-{cr_te_dists[-1]:.3f}]" if cr_te_dists else "Isolated"
    
    # 2. Point Group
    try:
        pg = SpacegroupAnalyzer(defect_sc, symprec=0.1).get_point_group_symbol()
    except:
        pg = "Unknown"
        
    # 3. Ewald Energy
    # Need to assign Cr oxidation state. Cr is usually +3 in these chalcogenides (isovalent with Sb+3).
    # Since Sb is +3 and Cr is +3, this is an isovalent substitution.
    # Charge neutral without compensation.
    defect_sc.add_oxidation_state_by_element({"Sb": 3, "Te": -2, "Cr": 3})
    
    try:
        es = EwaldSummation(defect_sc)
        ewald_energy = es.total_energy
    except Exception as e:
        ewald_energy = 0.0
        print(f"Ewald error: {e}")
        
    # Save
    filename = f"Cr_Sb_{wyckoff}_site_{idx}.xyz"
    filepath = os.path.join(output_dir, filename)
    
    atoms_out = AseAtomsAdaptor.get_atoms(defect_sc)
    atoms_out.info['total_charge'] = 0.0
    atoms_out.info['defect_type'] = "extrinsic_substitution"
    atoms_out.info['initial_magmom'] = 3.0 # For Cr3+ (d3 high spin? 3 muB)
    atoms_out.info['ewald_energy'] = ewald_energy
    atoms_out.info['point_group'] = pg
    
    # Set magnetic moments array for VASP
    mags = [0.0] * len(atoms_out)
    mags[idx] = 3.0
    atoms_out.set_initial_magnetic_moments(mags)
    
    ase_write(filepath, atoms_out, format="extxyz")
    
    generated_configs.append({
        "filename": filename,
        "type": "single",
        "site_index": idx,
        "wyckoff": wyckoff,
        "multiplicity": len(orbit['indices']),
        "point_group": pg,
        "nn_dist": min_dist,
        "nn_summary": nn_summary,
        "ewald_energy": ewald_energy
    })

# B. Multiple Substitution (x=0.05 approx)
# Target: ~3 Cr atoms on Sb sites (3/96 ~ 3.1%)
n_cr_target = 3
print(f"Generating multi-substitution configuration ({n_cr_target} Cr atoms)...")

# We want diverse spatial sampling (Maximin)
all_sb_indices = []
for o in sb_orbits:
    all_sb_indices.extend(o['indices'])

selected_indices = []
first_idx = random.choice(all_sb_indices)
selected_indices.append(first_idx)

for _ in range(n_cr_target - 1):
    best_idx = -1
    max_min_dist = -1.0
    
    for candidate in all_sb_indices:
        if candidate in selected_indices: continue
        
        pos_c = supercell[candidate].coords
        min_d = float('inf')
        for sel in selected_indices:
            pos_s = supercell[sel].coords
            d = np.linalg.norm(pos_c - pos_s)
            if d < min_d:
                min_d = d
        
        if min_d > max_min_dist:
            max_min_dist = min_d
            best_idx = candidate
            
    if best_idx != -1:
        selected_indices.append(best_idx)

defect_multi = supercell.copy()
for idx in selected_indices:
    defect_multi.replace(idx, "Cr")

# Properties
defect_multi.add_oxidation_state_by_element({"Sb": 3, "Te": -2, "Cr": 3})
try:
    ewald_multi = EwaldSummation(defect_multi).total_energy
except:
    ewald_multi = 0.0

filename_multi = f"Cr_Sb_Multi_x0.05.xyz"
filepath_multi = os.path.join(output_dir, filename_multi)

atoms_multi = AseAtomsAdaptor.get_atoms(defect_multi)
atoms_multi.info['total_charge'] = 0.0
atoms_multi.info['defect_type'] = "extrinsic_substitution_multi"
atoms_multi.info['concentration'] = 0.05
atoms_multi.info['ewald_energy'] = ewald_multi

mags_multi = [0.0] * len(atoms_multi)
for idx in selected_indices:
    mags_multi[idx] = 3.0
atoms_multi.set_initial_magnetic_moments(mags_multi)

ase_write(filepath_multi, atoms_multi, format="extxyz")

generated_configs.append({
    "filename": filename_multi,
    "type": "multi",
    "num_dopants": n_cr_target,
    "site_indices": selected_indices,
    "ewald_energy": ewald_multi,
    "note": "Maximin distribution"
})

# Output Report
# Sort single configs by Ewald Energy (Stability)
# (Though for isovalent, differences might be small/zero if symmetry equivalent)
single_configs = [c for c in generated_configs if c['type'] == 'single']
single_configs.sort(key=lambda x: x['ewald_energy'])

summary_report = {
    "summary": "Systematic Substitutional Cr->Sb Generation",
    "host": "Sb2Te3 4x4x1 Supercell",
    "dopant": "Cr (Substituting Sb)",
    "concentration_multi": "Approx x=0.05 (3 atoms / 96 sites)",
    "configurations": generated_configs
}

report_path = os.path.join(output_dir, "substitution_report.json")
with open(report_path, 'w') as f:
    json.dump(summary_report, f, indent=4)

print(f"Generation complete. Report: {report_path}")
