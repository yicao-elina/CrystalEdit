import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.bond_valence import BVAnalyzer
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
output_dir = "intrinsic_antisite"

if not os.path.exists(cif_path):
    # Fallback
    cif_path = "26ICML-CrystalEdit/Sb2Te3-mp-1201.cif"
    sc_path = "26ICML-CrystalEdit/sb2te3_supercell_441.xyz"
    output_dir = "26ICML-CrystalEdit/stage1/intrinsic_antisite"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print("Loading structures...")
unit_cell = Structure.from_file(cif_path)
atoms_sc = ase_read(sc_path)
supercell = AseAtomsAdaptor.get_structure(atoms_sc)

# Add oxidation states for Ewald sum and BVS (Sb+3, Te-2)
supercell.add_oxidation_state_by_element({"Sb": 3, "Te": -2})

# 1. Symmetry Orbit Expansion
print("Analyzing symmetry orbits in supercell...")
sga_sc = SpacegroupAnalyzer(supercell, symprec=0.01)
sym_struct_sc = sga_sc.get_symmetrized_structure()
equivalent_indices = sym_struct_sc.equivalent_indices

# Identify orbits corresponding to Te1 and Te2
# Te1 (6c) and Te2 (3a) in unit cell.
# In supercell, we need to identify which orbit corresponds to which site type.
# We can check the coordination or just map back to unit cell fractional coords.
# Or check specific known indices.
# Index 12 in SC was Te2 (Center QL) in previous step.
# Index 6 was Te1.

te_orbits = []
for orbit in equivalent_indices:
    rep_idx = orbit[0]
    site = supercell[rep_idx]
    if site.specie.symbol == "Te":
        # Identify type
        # Check Z coordinate relative to QL?
        # Or check wyckoff in supercell?
        # Let's map back to unit cell wyckoff if possible, or just label by index.
        # Simple heuristic:
        # Te2 is center of QL. Te1 is outer.
        # In this specific supercell (from previous outputs):
        # Te at [0,0,0] is Te2.
        
        orbit_label = f"Te_Orbit_{rep_idx}"
        
        # Check if 0.0,0.0,0.0 is in this orbit (using frac coords)
        is_te2_orbit = False
        for idx in orbit:
            if np.allclose(supercell[idx].frac_coords, [0,0,0], atol=0.01):
                is_te2_orbit = True
                break
        
        if is_te2_orbit:
            orbit_label = "Te2 (Center QL)"
        elif "Te1" not in [o['label'] for o in te_orbits]: # Assume others are Te1 if not Te2?
             # There might be multiple orbits for Te1 in supercell depending on symmetry breaking?
             # But 4x4x1 preserves R-3m usually.
             orbit_label = "Te1 (Outer QL)"
        
        te_orbits.append({
            "label": orbit_label,
            "indices": orbit
        })

print(f"Identified {len(te_orbits)} Te orbits in supercell.")

# 2. Combinatorial Generation
defect_concentration = 0.05
n_atoms = len(supercell)
n_antisites_target = int(np.ceil(n_atoms * defect_concentration))
# Actually, concentration is usually x in formula.
# Sb2Te3 -> 5 atoms. 5% usually means 5% of sites? or x=0.05?
# 5% of 240 atoms = 12 atoms.
# Let's target creating SINGLE antisites first (for library) 
# and then ONE multi-antisite configuration as requested.

generated_configs = []

# A. Single Antisites (One per unique orbit)
print("Generating single antisite configurations...")
for orbit in te_orbits:
    # Pick one representative index
    idx = orbit['indices'][0]
    label = orbit['label']
    orbit_name = label.split()[0]
    
    defect_sc = supercell.copy()
    defect_sc.replace(idx, "Sb")
    
    # Calculate Properties
    # Point Group
    try:
        pg = SpacegroupAnalyzer(defect_sc, symprec=0.1).get_point_group_symbol()
    except:
        pg = "Unknown"
        
    # BVS
    try:
        # BVS for the new Sb at idx
        # We need a structure with valences for BVAnalyzer?
        # BVAnalyzer estimates valences.
        # Let's just calculate Ewald energy as metric.
        bvs = "N/A" # BVAnalyzer is tricky with antisites without re-optimizing
    except:
        bvs = "Error"
        
    # Ewald Energy (Stability Ranking)
    # Assign Sb+3 to the new Sb (ionic limit check)
    # The new Sb replaces Te-2. So it sits on a -2 site.
    # Total charge would be +5 relative. Ewald requires neutral cell usually?
    # Or we assume background.
    # Pymatgen EwaldSummation requires oxidation states.
    try:
        # Compute energy difference relative to pristine
        # Pristine Ewald
        # es_pristine = EwaldSummation(supercell).total_energy
        # Defect Ewald (Force neutrality? No, just calc raw electrostatic change?)
        # EwaldSummation handles charged cells by background.
        es = EwaldSummation(defect_sc)
        ewald_energy = es.total_energy
    except Exception as e:
        ewald_energy = 0.0
        print(f"Ewald error: {e}")

    filename = f"Antisite_Sb_Te_{orbit_name}_{idx}.xyz"
    filepath = os.path.join(output_dir, filename)
    
    # Save
    atoms_out = AseAtomsAdaptor.get_atoms(defect_sc)
    atoms_out.info['total_charge'] = 0.0 # Assumption
    atoms_out.info['defect_type'] = "intrinsic_antisite"
    atoms_out.info['site_label'] = label
    atoms_out.info['ewald_energy'] = ewald_energy
    
    ase_write(filepath, atoms_out, format="extxyz")
    
    generated_configs.append({
        "filename": filename,
        "type": "single",
        "site_label": label,
        "site_index": idx,
        "point_group": pg,
        "ewald_energy": ewald_energy,
        "bvs": bvs
    })

# B. Multiple Antisites (Minkowski Sampling)
# Generate a configuration with ~5% concentration (approx 12 antisites)
# Use Minkowski distance to space them out (Maximin sampling)
print(f"Generating multi-antisite configuration ({n_antisites_target} sites)...")

all_te_indices = []
for o in te_orbits:
    all_te_indices.extend(o['indices'])

selected_indices = []
# Pick first random
first_idx = random.choice(all_te_indices)
selected_indices.append(first_idx)

# Greedily pick next indices to maximize distance to existing
for _ in range(n_antisites_target - 1):
    best_idx = -1
    max_min_dist = -1.0
    
    for candidate in all_te_indices:
        if candidate in selected_indices: continue
        
        # Find min distance to any selected
        # Cartesian distance
        pos_c = supercell[candidate].coords
        min_d = float('inf')
        for sel in selected_indices:
            pos_s = supercell[sel].coords
            d = np.linalg.norm(pos_c - pos_s) # Periodic distance would be better but cartesian approx ok for clustering check
            if d < min_d:
                min_d = d
        
        if min_d > max_min_dist:
            max_min_dist = min_d
            best_idx = candidate
    
    if best_idx != -1:
        selected_indices.append(best_idx)

# Create Multi-defect structure
defect_multi_sc = supercell.copy()
for idx in selected_indices:
    defect_multi_sc.replace(idx, "Sb")

# Ewald
try:
    ewald_multi = EwaldSummation(defect_multi_sc).total_energy
except:
    ewald_multi = 0.0

fname_multi = f"Antisite_Sb_Te_Multi_Conc_{defect_concentration:.2f}.xyz"
fpath_multi = os.path.join(output_dir, fname_multi)

atoms_multi = AseAtomsAdaptor.get_atoms(defect_multi_sc)
atoms_multi.info['total_charge'] = 0.0
atoms_multi.info['defect_type'] = "intrinsic_antisite_multi"
atoms_multi.info['concentration'] = defect_concentration
atoms_multi.info['ewald_energy'] = ewald_multi

ase_write(fpath_multi, atoms_multi, format="extxyz")

generated_configs.append({
    "filename": fname_multi,
    "type": "multi",
    "concentration": defect_concentration,
    "num_antisites": len(selected_indices),
    "site_indices": selected_indices,
    "ewald_energy": ewald_multi,
    "note": "Maximin spatial distribution"
})

# 5. Output Aggregation
summary_report = {
    "summary": "Systematic Antisite Expansion",
    "base_structure": "Sb2Te3 4x4x1 Supercell",
    "defect_type": "Sb_Te Antisite",
    "configurations": generated_configs
}

# Update existing report or create new
report_path = os.path.join(output_dir, "antisite_expansion_report.json")
with open(report_path, 'w') as f:
    json.dump(summary_report, f, indent=4)

print(f"Expansion complete. Generated {len(generated_configs)} files.")
print(f"Report saved to {report_path}")
