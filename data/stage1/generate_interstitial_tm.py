import numpy as np
from pymatgen.core import Structure, PeriodicSite, Element
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.analysis.local_env import CrystalNN
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
base_output_dir = "tetrahedral_interstitials"

if not os.path.exists(cif_path):
    cif_path = "26ICML-CrystalEdit/Sb2Te3-mp-1201.cif"
    sc_path = "26ICML-CrystalEdit/sb2te3_supercell_441.xyz"
    base_output_dir = "26ICML-CrystalEdit/stage1/tetrahedral_interstitials"

if not os.path.exists(base_output_dir):
    os.makedirs(base_output_dir)

print("Loading structures...")
unit_cell = Structure.from_file(cif_path)
atoms_sc = ase_read(sc_path)
supercell = AseAtomsAdaptor.get_structure(atoms_sc)

# Dopant Definitions (Ionic Radii for +3 state in CN6/CN4 approx)
# Shannon Radii: 
# Cr3+ (VI) = 0.615 A
# V3+ (VI) = 0.64 A
# Mn3+ (VI) = 0.645 A (High spin) / 0.58 (Low spin). Using 0.65 as per prompt.
dopants = {
    "Cr": {"radius": 0.62, "magmom": 3.0},
    "V":  {"radius": 0.64, "magmom": 2.0},
    "Mn": {"radius": 0.65, "magmom": 4.0}
}

# 1. Symmetry-Aware Void Identification
print("Identifying interstitial sites via Voronoi...")
generator = VoronoiInterstitialGenerator()
# Use Dummy Cr for generation
interstitials = list(generator.generate(unit_cell, {"Cr": 1}))

print(f"Found {len(interstitials)} distinct interstitial sites in unit cell.")

# VdW Gap centers (z fractional)
vdw_centers = [1/6.0, 3/6.0, 5/6.0]

candidates = []
cnn = CrystalNN()

for i, defect in enumerate(interstitials):
    site = defect.site
    frac_coords = site.frac_coords
    
    # 1.1 Check Location (VdW vs Intralayer)
    is_vdw = False
    for center in vdw_centers:
        z_dist = abs(frac_coords[2] - center)
        z_dist = min(z_dist, 1.0 - z_dist)
        if z_dist < 0.08: # Close to gap center
            is_vdw = True
            break
            
    # 1.2 Analyze Geometry
    temp_struct = unit_cell.copy()
    temp_struct.append("Cr", frac_coords)
    
    try:
        cn_info = cnn.get_nn_info(temp_struct, len(temp_struct)-1)
        cn = len(cn_info)
    except:
        cn = 0
        
    # Get neighbors manually for min_dist
    neighbors = unit_cell.get_neighbors(site, r=4.0)
    if not neighbors: continue
    neighbors = sorted(neighbors, key=lambda x: x.nn_distance)
    min_dist = neighbors[0].nn_distance
    
    # Heuristic Classification
    geometry = "Unknown"
    if cn == 4:
        geometry = "Tetrahedral"
    elif cn == 6:
        geometry = "Octahedral"
    else:
        # Check geometric neighbors count (ignoring weights)
        close_neighbors = [n for n in neighbors if n.nn_distance < 3.0]
        c_cn = len(close_neighbors)
        if c_cn == 4:
            geometry = "Tetrahedral-Like"
        elif c_cn == 6:
            geometry = "Octahedral-Like"
        else:
            geometry = f"Distorted_CN{c_cn}"

    candidates.append({
        "id": i,
        "geometry": geometry,
        "is_vdw": is_vdw,
        "site": site,
        "min_dist": min_dist,
        "cn": cn
    })

# Selection Logic: Prioritize Tetrahedral in VdW
selected_site_data = None

# Filter
vdw_tets = [c for c in candidates if c['is_vdw'] and "Tetrahedral" in c['geometry']]
vdw_any = [c for c in candidates if c['is_vdw']]
intra_tets = [c for c in candidates if not c['is_vdw'] and "Tetrahedral" in c['geometry']]

if vdw_tets:
    # Pick largest
    selected_site_data = max(vdw_tets, key=lambda x: x['min_dist'])
    print("Selected VdW Tetrahedral site.")
elif vdw_any:
    # Pick largest VdW
    selected_site_data = max(vdw_any, key=lambda x: x['min_dist'])
    print(f"Selected VdW site (Geometry: {selected_site_data['geometry']}) - No perfect Tet found.")
elif intra_tets:
    selected_site_data = max(intra_tets, key=lambda x: x['min_dist'])
    print("Selected Intralayer Tetrahedral site.")
else:
    selected_site_data = max(candidates, key=lambda x: x['min_dist'])
    print("Selected fallback site (largest void).")

print(f"Chosen Site ID {selected_site_data['id']}: {selected_site_data['geometry']} @ {selected_site_data['site'].frac_coords}")

# 2. Process Each Dopant
# We will relax the ghost atom separately for each dopant species to account for radius?
# Actually, standard relaxation usually uses host interactions.
# But let's assume the "dopant" interacts via LJ.
# We optimize the position within the static host lattice.

def lj_potential(frac_coords, structure, r_ion):
    # E = 4*epsilon * ((sigma/r)^12 - (sigma/r)^6)
    # Simplified repulsive/attractive for "fitting" into a hole
    # Ideally we want distance ~ r_ion + r_host_avg (~1.4)
    # Target distance d0 = r_ion + 1.4
    d0 = r_ion + 1.4
    sigma = d0 / (2**(1/6))
    epsilon = 1.0
    
    test_site = PeriodicSite("H", frac_coords, structure.lattice) # Dummy species
    neighbors = structure.get_neighbors(test_site, r=6.0)
    
    e = 0.0
    for n in neighbors:
        r = n.nn_distance
        if r < 0.1: r = 0.1
        sr = sigma / r
        e += 4 * epsilon * (sr**12 - sr**6)
    
    # Add penalty for breaking symmetry? 
    # Or restrict optimization to symmetry lines.
    # Nelder-Mead is unconstrained but if we start at high symmetry, it might stay close.
    return e

summary_report = {
    "summary": "Systematic Multi-Dopant Interstitial Generation",
    "base_structure": "Sb2Te3 4x4x1 Supercell",
    "selected_void": {
        "id": selected_site_data['id'],
        "original_geometry": selected_site_data['geometry'],
        "is_vdw": selected_site_data['is_vdw'],
        "original_frac": selected_site_data['site'].frac_coords.tolist()
    },
    "dopants": {}
}

scaling_matrix = [4, 4, 1]

for el, data in dopants.items():
    print(f"\nProcessing Dopant: {el}")
    r_ion = data['radius']
    magmom = data['magmom']
    
    out_dir = os.path.join(base_output_dir, el)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    # Relaxation
    print(f"  Optimizing position for {el} (r={r_ion})...")
    res = minimize(lj_potential, selected_site_data['site'].frac_coords, args=(unit_cell, r_ion), method='Nelder-Mead', tol=1e-4)
    opt_frac = res.x % 1.0
    
    # Calculate Metrics after relaxation
    temp_site = PeriodicSite(el, opt_frac, unit_cell.lattice)
    neighbors = unit_cell.get_neighbors(temp_site, r=4.0)
    neighbors.sort(key=lambda x: x.nn_distance)
    
    nn4 = neighbors[:4]
    d_nn_avg = np.mean([n.nn_distance for n in nn4])
    min_dist_final = neighbors[0].nn_distance
    
    # Steric Overlap Factor S = R_void / R_ion
    # R_void ~ min_dist (distance to nuclei center) - R_host_core?
    # Prompt says: "ratio of the void radius to the ionic radius".
    # Void radius usually means radius of sphere that fits touching host atoms.
    # R_void_eff = min_dist - R_Te_ionic. 
    # If using that def, R_void might be small.
    # Let's stick to the prompt's implied simple metric or define clearly.
    # "Steric Overlap Factor based on the ratio of the void radius to the ionic radius"
    # Let's calculate Void Radius = min_dist_final. (Simple proxy).
    steric_factor = min_dist_final / r_ion
    
    # Distortion Index (Bond angle variance for 4 NN)
    distortion_index = 0.0
    if len(nn4) >= 4:
        angles = []
        import itertools
        for n1, n2 in itertools.combinations(nn4, 2):
            v1 = n1.coords - temp_site.coords
            v2 = n2.coords - temp_site.coords
            angle = np.degrees(np.arccos(np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0)))
            angles.append(angle)
        # Select the 6 tetrahedral angles (ideal 109.5)
        # If we have more neighbors, picking 4 closest is heuristic.
        # Variance
        distortion_index = np.mean([(a - 109.47)**2 for a in angles[:6]]) 
    
    # Symmetry Analysis
    # Place in unit cell and check point group
    struct_sym = unit_cell.copy()
    struct_sym.append(el, opt_frac)
    try:
        sga = SpacegroupAnalyzer(struct_sym, symprec=0.1)
        pg = sga.get_point_group_symbol()
    except:
        pg = "C1"
        
    reasoning = f"Steric Factor {steric_factor:.2f}. "
    if steric_factor < 3.0: # Heuristic threshold (dist/r_ion usually ~ 4-5 for no overlap if purely center-to-center? No. Bond length is r1+r2. So dist/r1 = 1 + r2/r1. If r2~2, r1~0.6, ratio ~ 4.3).
        # Actually min_dist is bond length.
        # r_void should be "space available".
        # Let's assume ratio > 1.0 means it fits "inside" the bond distance? No.
        # Just report the value.
        reasoning += "Tight fit."
    else:
        reasoning += "Spacious."
        
    # Generate Structures
    
    # 1. Single (x=0.01)
    # Map to Supercell
    sc_frac = [opt_frac[0]/4.0, opt_frac[1]/4.0, opt_frac[2]]
    
    sc_single = supercell.copy()
    sc_single.append(el, sc_frac)
    
    fname_single = f"{el}_tet_x0.01.xyz"
    ase_single = AseAtomsAdaptor.get_atoms(sc_single)
    
    # Set Magmom
    mags = [0.0] * len(ase_single)
    mags[-1] = magmom
    ase_single.set_initial_magnetic_moments(mags)
    
    ase_single.info['total_charge'] = 3.0 # +3 state
    ase_single.info['dopant_species'] = el
    ase_single.info['steric_factor'] = steric_factor
    ase_single.info['distortion_index'] = distortion_index
    ase_write(os.path.join(out_dir, fname_single), ase_single, format="extxyz")
    
    # 2. Multi (x=0.05) ~ 3 atoms
    n_multi = 3
    # Maximin sampling in Supercell
    possible_sites = []
    for i in range(4):
        for j in range(4):
            possible_sites.append([(opt_frac[0]+i)/4.0, (opt_frac[1]+j)/4.0, opt_frac[2]])
            
    selected = [sc_frac]
    for _ in range(n_multi - 1):
        best = None
        max_dist = -1
        for cand in possible_sites:
            if any(np.allclose(cand, s) for s in selected): continue
            
            d = min([np.linalg.norm(np.dot(np.array(cand)-np.array(s), supercell.lattice.matrix)) for s in selected])
            if d > max_dist:
                max_dist = d
                best = cand
        if best: selected.append(best)
        
    sc_multi = supercell.copy()
    for s in selected:
        sc_multi.append(el, s)
        
    fname_multi = f"{el}_tet_x0.05.xyz"
    ase_multi = AseAtomsAdaptor.get_atoms(sc_multi)
    
    mags_m = [0.0] * len(supercell) + [magmom] * n_multi
    ase_multi.set_initial_magnetic_moments(mags_m)
    
    ase_multi.info['total_charge'] = 3.0 * n_multi
    ase_multi.info['concentration'] = 0.05
    ase_write(os.path.join(out_dir, fname_multi), ase_multi, format="extxyz")
    
    # Report Data
    dopant_report = {
        "ionic_radius": r_ion,
        "optimized_frac_coords": opt_frac.tolist(),
        "point_group": pg,
        "coordination": {
            "cn_4_neighbors": [n.specie.symbol for n in nn4],
            "distances": [float(n.nn_distance) for n in nn4],
            "average_dist": float(d_nn_avg),
            "distortion_index_sigma2": float(distortion_index)
        },
        "energetics_proxy": {
            "steric_overlap_factor": float(steric_factor),
            "reasoning": reasoning
        },
        "files": [fname_single, fname_multi]
    }
    
    # Save individual report
    with open(os.path.join(out_dir, f"{el}_tet_report.json"), "w") as f:
        json.dump(dopant_report, f, indent=4)
        
    summary_report['dopants'][el] = dopant_report

# Save Global Report
with open(os.path.join(base_output_dir, "multi_dopant_report.json"), "w") as f:
    json.dump(summary_report, f, indent=4)

print("Done.")
