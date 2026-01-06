import numpy as np
from pymatgen.core import Structure, PeriodicSite
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.analysis.local_env import CrystalNN, LocalStructOrderParams
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read as ase_read
from ase.io import write as ase_write
import json
import os
import random
from scipy.optimize import minimize

# Paths
cif_path = "../Sb2Te3-mp-1201.cif"
sc_path = "../sb2te3_supercell_441.xyz"
output_dir = "interstitial_hydrogen"

if not os.path.exists(cif_path):
    # Fallback
    cif_path = "26ICML-CrystalEdit/Sb2Te3-mp-1201.cif"
    sc_path = "26ICML-CrystalEdit/sb2te3_supercell_441.xyz"
    output_dir = "26ICML-CrystalEdit/stage1/interstitial_hydrogen"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print("Loading structures...")
unit_cell = Structure.from_file(cif_path)
atoms_sc = ase_read(sc_path)
supercell = AseAtomsAdaptor.get_structure(atoms_sc)

# 1. Symmetry-Aware Void Identification
print("Identifying interstitial sites via Voronoi...")
# Use unit cell for identification to find unique types, then map to supercell
generator = VoronoiInterstitialGenerator()
# Insert H
interstitials = list(generator.generate(unit_cell, {"H": 1}))

print(f"Found {len(interstitials)} distinct interstitial sites in unit cell.")

# Helper to identify VdW gap
# In Sb2Te3 (R-3m) with 3 formula units, gaps are at z ~ 1/6, 3/6, 5/6
vdw_centers = [1/6.0, 3/6.0, 5/6.0]
print(f"Assumed VdW gap centers at Z = {vdw_centers}")

candidates = []
cnn = CrystalNN()

for i, defect in enumerate(interstitials):
    site = defect.site
    # Analyze Coordination
    temp_struct = unit_cell.copy()
    temp_struct.append("H", site.frac_coords)
    
    # Nearest Neighbors
    # Get neighbors within 3.5 A
    neighbors = unit_cell.get_neighbors(site, r=3.5)
    neighbors = sorted(neighbors, key=lambda x: x.nn_distance)
    
    if not neighbors: continue
    
    min_dist = neighbors[0].nn_distance
    cn = len(neighbors) # Simple CN
    
    # CrystalNN for better CN
    try:
        cn_info = cnn.get_nn_info(temp_struct, len(temp_struct)-1)
        real_cn = len(cn_info)
    except:
        real_cn = 0
        
    # Fallback if CrystalNN fails to find neighbors for small H
    if real_cn == 0:
        real_cn = cn
    
    # Classification
    # Check if in VdW gap (Z coordinate)
    is_vdw = False
    for center in vdw_centers:
        z_dist = abs(site.frac_coords[2] - center)
        z_dist = min(z_dist, 1.0 - z_dist) # PBC
        if z_dist < 0.08: # Approx half width of gap
            is_vdw = True
            break
    
    geometry = "Unknown"
    if real_cn == 4:
        geometry = "Tetrahedral"
    elif real_cn == 6:
        geometry = "Octahedral"
    else:
        geometry = f"Distorted_CN{real_cn}"
        
    candidates.append({
        "id": i,
        "geometry": geometry,
        "is_vdw": is_vdw,
        "site": site,
        "min_dist": min_dist,
        "cn": real_cn
    })

# Selection Logic: Prioritize Tetrahedral in VdW
selected_candidate = None

# Filter for VdW Tetrahedral
vdw_tets = [c for c in candidates if c['is_vdw'] and "Tetrahedral" in c['geometry']]
vdw_any = [c for c in candidates if c['is_vdw']]
intra_tets = [c for c in candidates if not c['is_vdw'] and "Tetrahedral" in c['geometry']]

print(f"Sites found: {len(vdw_tets)} VdW-Tet, {len(vdw_any)} VdW-Any, {len(intra_tets)} Intra-Tet")

if vdw_tets:
    # Pick largest (max min_dist)
    selected_candidate = max(vdw_tets, key=lambda x: x['min_dist'])
    print("Selected VdW Tetrahedral site.")
elif vdw_any:
    # Fallback to any VdW
    # Prefer lower CN (closer to Tet)
    selected_candidate = min(vdw_any, key=lambda x: x['cn'])
    print(f"Selected VdW site (Geometry: {selected_candidate['geometry']}).")
elif intra_tets:
    selected_candidate = max(intra_tets, key=lambda x: x['min_dist'])
    print("Selected Intralayer Tetrahedral site (No VdW found).")
else:
    selected_candidate = candidates[0]
    print("Selected fallback site.")

print(f"Chosen Site ID {selected_candidate['id']}: {selected_candidate['geometry']} @ {selected_candidate['site'].frac_coords}")

# 2. Local Relaxation (Geometric)
# H is small. We place it. We can perform a simple "nudging" to maximize distance to neighbors.
# Define potential: U = sum (1/r^12) repulsive
def potential(frac_coords, structure):
    test_site = PeriodicSite("H", frac_coords, structure.lattice)
    neighbors = structure.get_neighbors(test_site, r=4.0)
    e = 0.0
    for n in neighbors:
        d = n.nn_distance
        if d < 0.1: d = 0.1
        e += (2.0 / d)**12 # Simple repulsion, scale 2.0A
    return e

# Optimize in Unit Cell
# Constraints: Keep symmetry? 
# If we start at Voronoi node (high symmetry typically), gradient might be zero.
# Let's try to minimize.
res = minimize(potential, selected_candidate['site'].frac_coords, args=(unit_cell), method='Nelder-Mead', tol=1e-3)
optimized_frac = res.x
optimized_frac = optimized_frac % 1.0

# Update selected site
old_site = selected_candidate['site']
selected_candidate['site'] = PeriodicSite("H", optimized_frac, unit_cell.lattice)
shift_cart = np.linalg.norm(selected_candidate['site'].coords - old_site.coords)
print(f"Optimized Position: {optimized_frac} (Cartesian Shift: {shift_cart:.4f} A)")

# 3. Map to Supercell and Generate Single Interstitial
# We need to pick ONE location in the supercell corresponding to this unit cell site.
# 4x4x1
target_sc_frac = optimized_frac.copy()
target_sc_frac[0] /= 4.0
target_sc_frac[1] /= 4.0
# Z is same

defect_sc = supercell.copy()
defect_sc.append("H", target_sc_frac)

# Metrics
# Tetrahedral Distortion Index (if CN=4)
# Calculate bond angle variance
neighbors_sc = defect_sc.get_neighbors(defect_sc[-1], r=3.5)
neighbors_sc = sorted(neighbors_sc, key=lambda x: x.nn_distance)
nn4 = neighbors_sc[:4]
distortion_index = "N/A"
if len(nn4) == 4:
    angles = []
    import itertools
    for a, b in itertools.combinations(nn4, 2):
        # Calculate angle H-a-b? No, angle a-H-b
        # H is at center
        v1 = a.coords - defect_sc[-1].coords
        v2 = b.coords - defect_sc[-1].coords
        angle = np.degrees(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
        angles.append(angle)
    avg_angle = np.mean(angles)
    variance = np.sum([(theta - 109.5)**2 for theta in angles]) / 6.0 # Ideal is 109.5
    distortion_index = f"{variance:.2f}"

# Steric Overlap Factor
# R_void ~ min_dist. R_host_cov ~ 1.4 (Sb/Te avg)
steric_factor = selected_candidate['min_dist'] / 1.4

# Save Single
filename_single = f"H_Interstitial_{selected_candidate['geometry']}_{selected_candidate['id']}.xyz"
filepath_single = os.path.join(output_dir, filename_single)

atoms_single = AseAtomsAdaptor.get_atoms(defect_sc)
atoms_single.info['total_charge'] = 1.0 # H+
atoms_single.info['defect_type'] = "interstitial_H"
atoms_single.info['distortion_index'] = distortion_index
atoms_single.info['steric_factor'] = steric_factor

ase_write(filepath_single, atoms_single, format="extxyz")

# 4. Multi-H Generation (x=0.05)
# Target ~ 2-3 atoms. Let's do 3.
n_h = 3
print(f"Generating multi-H configuration ({n_h} atoms)...")

# Generate all possible sites in supercell
possible_sc_sites = []
for i in range(4):
    for j in range(4):
        # Apply shift
        f = optimized_frac.copy()
        sc_f = [(f[0]+i)/4.0, (f[1]+j)/4.0, f[2]]
        possible_sc_sites.append(sc_f)

# Maximin Sampling
selected_sites = []
first_site = random.choice(possible_sc_sites)
selected_sites.append(first_site)

for _ in range(n_h - 1):
    best_site = None
    max_min_dist = -1.0
    
    for cand in possible_sc_sites:
        # Skip if already close to selected (identity check)
        if any(np.allclose(cand, s, atol=0.01) for s in selected_sites): continue
        
        # Calc distance to existing set
        cand_cart = np.dot(cand, supercell.lattice.matrix)
        min_d = float('inf')
        for sel in selected_sites:
            sel_cart = np.dot(sel, supercell.lattice.matrix)
            d = np.linalg.norm(cand_cart - sel_cart)
            if d < min_d: min_d = d
            
        if min_d > max_min_dist:
            max_min_dist = min_d
            best_site = cand
    
    if best_site is not None:
        selected_sites.append(best_site)

defect_multi = supercell.copy()
for s in selected_sites:
    defect_multi.append("H", s)

filename_multi = f"H_Interstitial_Multi_x0.05.xyz"
filepath_multi = os.path.join(output_dir, filename_multi)

atoms_multi = AseAtomsAdaptor.get_atoms(defect_multi)
atoms_multi.info['total_charge'] = float(n_h) # +1 per H
atoms_multi.info['defect_type'] = "interstitial_H_multi"
atoms_multi.info['concentration'] = 0.05

ase_write(filepath_multi, atoms_multi, format="extxyz")

# Report
report = {
    "summary": "Systematic Interstitial Hydrogen Generation",
    "site_selection": {
        "geometry": selected_candidate['geometry'],
        "is_vdw": bool(selected_candidate['is_vdw']),
        "min_dist_to_host": float(selected_candidate['min_dist']),
        "optimization_shift": float(np.linalg.norm(res.x - selected_candidate['site'].frac_coords))
    },
    "metrics": {
        "tetrahedral_distortion_variance": distortion_index,
        "steric_overlap_factor": float(steric_factor),
        "physical_reasoning": f"This site is {'meta-stable' if steric_factor > 0.8 else 'likely unstable'} because of the steric overlap factor ({steric_factor:.2f}). VdW placement reduces strain."
    },
    "configurations": [
        {"filename": filename_single, "type": "single", "charge": 1.0},
        {"filename": filename_multi, "type": "multi", "charge": float(n_h), "note": "Maximin sampling"}
    ]
}

report_path = os.path.join(output_dir, "interstitial_report.json")
with open(report_path, 'w') as f:
    json.dump(report, f, indent=4)

print(f"Completed. Report: {report_path}")
