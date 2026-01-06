import numpy as np
import csv
from pymatgen.core import Structure, PeriodicSite
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read as ase_read
from ase.io import write as ase_write
import json
import os
import itertools
from scipy.spatial import ConvexHull

# Definitions
dopants = {
    "Cr": {"radius": 0.62, "magmom": 3.0},
    "V":  {"radius": 0.64, "magmom": 2.0},
    "Mn": {"radius": 0.65, "magmom": 4.0}
}

# Paths
cif_path = "26ICML-CrystalEdit/Sb2Te3-mp-1201.cif"
sc_path = "26ICML-CrystalEdit/sb2te3_supercell_441.xyz"
base_output_dir = "26ICML-CrystalEdit/stage1/octahedral_interstitials"

if not os.path.exists(cif_path):
    # Fallback for relative paths
    cif_path = "../Sb2Te3-mp-1201.cif"
    sc_path = "../sb2te3_supercell_441.xyz"
    base_output_dir = "octahedral_interstitials"

if not os.path.exists(base_output_dir):
    os.makedirs(base_output_dir)

print("Loading structures...")
unit_cell = Structure.from_file(cif_path)
atoms_sc = ase_read(sc_path)
supercell = AseAtomsAdaptor.get_structure(atoms_sc)
# Add oxidation states for Ewald (Host)
supercell.add_oxidation_state_by_element({"Sb": 3, "Te": -2})

# 1. Identify Octahedral Voids
print("Identifying octahedral voids...")
generator = VoronoiInterstitialGenerator()
interstitials = list(generator.generate(unit_cell, {"Cr": 1}))

octahedral_candidates = []
vdw_centers = [1/6.0, 3/6.0, 5/6.0]

for i, defect in enumerate(interstitials):
    site = defect.site
    # Analyze Coordination
    temp_struct = unit_cell.copy()
    temp_struct.append("Cr", site.frac_coords)
    
    # Get neighbors
    neighbors = unit_cell.get_neighbors(site, r=4.0)
    neighbors.sort(key=lambda x: x.nn_distance)
    
    # Check location first
    is_vdw = False
    for center in vdw_centers:
        z_dist = abs(site.frac_coords[2] - center)
        z_dist = min(z_dist, 1.0 - z_dist)
        if z_dist < 0.08:
            is_vdw = True
            break

    # Check for CN=6
    # Heuristic: Look for 6 neighbors within a reasonable shell
    # Relaxed cutoff for VdW sites which might be larger
    close_neighbors = [n for n in neighbors if n.nn_distance < 3.6]
    
    is_octahedral = False
    if len(close_neighbors) >= 6:
        # Check if the gap between 6th and 7th neighbor is significant
        # or if the 6th neighbor is within valid bonding range
        d6 = neighbors[5].nn_distance
        d7 = neighbors[6].nn_distance if len(neighbors) > 6 else 10.0
        
        # If 6 neighbors are clustered and 7th is far
        if d7 - d6 > 0.2:
            is_octahedral = True
        # Or if we just have 6 close ones
        elif len(close_neighbors) == 6:
            is_octahedral = True
            
    # Always keep VdW candidates for inspection if CN is high
    if not is_octahedral:
        # If in VdW gap and CN>=6, considering it a candidate to relax
        if is_vdw and len(close_neighbors) >= 6:
            is_octahedral = True

    if is_octahedral:
        octahedral_candidates.append({
            "id": i,
            "is_vdw": is_vdw,
            "site": site,
            "neighbors": neighbors[:6]
        })

print(f"Found {len(octahedral_candidates)} octahedral candidates.")

# Selection: VdW Octahedral
selected_void = None
vdw_octs = [c for c in octahedral_candidates if c['is_vdw']]
intra_octs = [c for c in octahedral_candidates if not c['is_vdw']]

if vdw_octs:
    # Pick the one with most regular geometry (lowest variance in bond lengths)
    def variance(c):
        d = [n.nn_distance for n in c['neighbors']]
        return np.var(d)
    selected_void = min(vdw_octs, key=variance)
    print("Selected VdW Octahedral site.")
elif intra_octs:
    selected_void = intra_octs[0]
    print("Selected Intralayer Octahedral site (No VdW found).")
else:
    print("No octahedral sites found. Using fallback logic (Voronoi node).")
    # Fallback: Just pick the largest VdW void
    vdw_all = [i for i, d in enumerate(interstitials) if any(abs(d.site.frac_coords[2]-c)<0.08 for c in vdw_centers)]
    if vdw_all:
        site = interstitials[vdw_all[0]].site
        neighbors = unit_cell.get_neighbors(site, r=4.0)
        selected_void = {"site": site, "is_vdw": True, "id": vdw_all[0], "neighbors": neighbors[:6]}
    else:
        site = interstitials[0].site
        neighbors = unit_cell.get_neighbors(site, r=4.0)
        selected_void = {"site": site, "is_vdw": False, "id": 0, "neighbors": neighbors[:6]}

# 2. Geometry Metrics Functions
def calculate_metrics(site_coords, neighbors):
    # Neighbors is list of PeriodicSites
    coords_site = site_coords
    coords_nn = [n.coords for n in neighbors]
    
    # Bond lengths
    dists = [np.linalg.norm(c - coords_site) for c in coords_nn]
    avg_d = np.mean(dists)
    
    # Volume of coordination polyhedron
    try:
        hull = ConvexHull(coords_nn)
        vol = hull.volume
    except:
        vol = 0.0
        
    # Quadratic Elongation lambda
    # l_0 for ideal octahedron of same volume
    # V = 4/3 * l0^3 (where l0 is center-to-vertex distance)
    # l0 = (3V/4)^(1/3)
    if vol > 0.1:
        l0 = (3 * vol / 4.0)**(1/3)
        lambda_val = np.mean([(d/l0)**2 for d in dists])
    else:
        lambda_val = 1.0
        
    # Bond Angle Variance sigma^2
    angles = []
    for i in range(len(neighbors)):
        for j in range(i+1, len(neighbors)):
            v1 = coords_nn[i] - coords_site
            v2 = coords_nn[j] - coords_site
            angle = np.degrees(np.arccos(np.clip(np.dot(v1, v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)), -1.0, 1.0)))
            angles.append(angle)
    
    # For octahedron, ideal angles are 90 (12 angles) and 180 (3 angles).
    # Standard definition usually considers the 12 cis angles (90 deg).
    # Identify cis angles (approx 90)
    cis_angles = [a for a in angles if 70 < a < 110]
    if cis_angles:
        sigma2 = np.var([a - 90 for a in cis_angles]) # Mean square deviation from 90
        # Definition: sum(theta_i - 90)^2 / 11
        sigma2 = np.sum([(a-90)**2 for a in cis_angles]) / (len(cis_angles)-1) if len(cis_angles) > 1 else 0.0
    else:
        sigma2 = 0.0
        
    return avg_d, lambda_val, sigma2, vol

# 3. Processing
summary_rows = []

# Map selected unit cell void to supercell
uc_frac = selected_void['site'].frac_coords
sc_target_frac = [uc_frac[0]/4.0, uc_frac[1]/4.0, uc_frac[2]]

# Identify adjacent site for clustering
# Find another equivalent site in supercell closest to sc_target_frac
possible_sc_sites = []
for i in range(4):
    for j in range(4):
        # z is same for 4x4x1
        f = [(uc_frac[0]+i)/4.0, (uc_frac[1]+j)/4.0, uc_frac[2]]
        possible_sc_sites.append(f)

# Convert to cartesian to find closest
site0_cart = np.dot(sc_target_frac, supercell.lattice.matrix)
dists_to_0 = []
for f in possible_sc_sites:
    cart = np.dot(f, supercell.lattice.matrix)
    d = np.linalg.norm(cart - site0_cart)
    if d > 0.1: # Exclude self
        dists_to_0.append((d, f))
dists_to_0.sort(key=lambda x: x[0])
closest_site_frac = dists_to_0[0][1]

for el, props in dopants.items():
    print(f"\nGenerating for {el}...")
    r_ion = props['radius']
    magmom = props['magmom']
    
    out_dir = os.path.join(base_output_dir, el)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    
    # --- Single Dilute ---
    sc_dilute = supercell.copy()
    sc_dilute.append(el, sc_target_frac)
    # Add oxidation state for Ewald (Must provide all)
    sc_dilute.add_oxidation_state_by_element({"Sb": 3, "Te": -2, el: 3})
    
    # Ewald
    try:
        ewald_dilute = EwaldSummation(sc_dilute).total_energy
    except:
        ewald_dilute = 0.0
        
    # Metrics
    temp_site_obj = sc_dilute[-1]
    nn = sc_dilute.get_neighbors(temp_site_obj, r=4.0)
    nn.sort(key=lambda x: x.nn_distance)
    nn6 = nn[:6]
    
    avg_d, lam, sig2, vol = calculate_metrics(temp_site_obj.coords, nn6)
    
    # Symmetry
    try:
        pg = SpacegroupAnalyzer(sc_dilute, symprec=0.1).get_point_group_symbol()
    except:
        pg = "C1"
        
    # Delta V (Proxy: Vol of poly vs Ideal Vol for Radius?) 
    # Or just use Vol from ConvexHull
    
    # Save Dilute
    atoms_dilute = AseAtomsAdaptor.get_atoms(sc_dilute)
    mags = [0.0]*len(supercell) + [magmom]
    atoms_dilute.set_initial_magnetic_moments(mags)
    atoms_dilute.info.update({
        "dopant_species": el,
        "site_symmetry": pg,
        "initial_magmom": magmom,
        "config_type": "dilute",
        "ewald_sum": ewald_dilute
    })
    ase_write(os.path.join(out_dir, f"{el}_oct_x0.01_dilute.xyz"), atoms_dilute, format="extxyz")
    
    summary_rows.append({
        "dopant": el,
        "concentration": 0.01,
        "config_type": "dilute",
        "avg_d_nn": f"{avg_d:.3f}",
        "lambda": f"{lam:.4f}",
        "sigma2": f"{sig2:.2f}",
        "delta_V_poly": f"{vol:.2f}",
        "ewald_sum": f"{ewald_dilute:.2f}",
        "symmetry": pg
    })
    
    # --- Clustered (2 atoms) ---
    sc_cluster = supercell.copy()
    sc_cluster.append(el, sc_target_frac)
    sc_cluster.append(el, closest_site_frac)
    sc_cluster.add_oxidation_state_by_element({"Sb": 3, "Te": -2, el: 3})
    
    try:
        ewald_cluster = EwaldSummation(sc_cluster).total_energy
    except:
        ewald_cluster = 0.0
        
    atoms_cluster = AseAtomsAdaptor.get_atoms(sc_cluster)
    mags_c = [0.0]*len(supercell) + [magmom, magmom]
    atoms_cluster.set_initial_magnetic_moments(mags_c)
    atoms_cluster.info.update({
        "dopant_species": el,
        "config_type": "clustered",
        "ewald_sum": ewald_cluster
    })
    ase_write(os.path.join(out_dir, f"{el}_oct_x0.01_clustered.xyz"), atoms_cluster, format="extxyz")
    
    summary_rows.append({
        "dopant": el,
        "concentration": 0.01,
        "config_type": "clustered",
        "avg_d_nn": "N/A",
        "lambda": "N/A",
        "sigma2": "N/A",
        "delta_V_poly": "N/A",
        "ewald_sum": f"{ewald_cluster:.2f}",
        "symmetry": "N/A"
    })
    
    # --- Multi x=0.05 (3 atoms) ---
    # Maximin
    selected = [sc_target_frac]
    for _ in range(2): # Total 3
        best = None
        max_dist = -1
        for cand in possible_sc_sites:
            if any(np.allclose(cand, s) for s in selected): continue
            d = min([np.linalg.norm(np.dot(np.array(cand)-np.array(s), supercell.lattice.matrix)) for s in selected])
            if d > max_dist:
                max_dist = d
                best = cand
        if best: selected.append(best)
        
    sc_multi = supercell.copy()
    for s in selected:
        sc_multi.append(el, s)
    sc_multi.add_oxidation_state_by_element({"Sb": 3, "Te": -2, el: 3})
    
    try:
        ewald_multi = EwaldSummation(sc_multi).total_energy
    except:
        ewald_multi = 0.0
        
    atoms_multi = AseAtomsAdaptor.get_atoms(sc_multi)
    mags_m = [0.0]*len(supercell) + [magmom]*3
    atoms_multi.set_initial_magnetic_moments(mags_m)
    atoms_multi.info.update({
        "dopant_species": el,
        "config_type": "multi_x0.05",
        "ewald_sum": ewald_multi
    })
    ase_write(os.path.join(out_dir, f"{el}_oct_x0.05.xyz"), atoms_multi, format="extxyz")
    
    summary_rows.append({
        "dopant": el,
        "concentration": 0.05,
        "config_type": "multi",
        "avg_d_nn": "N/A",
        "lambda": "N/A",
        "sigma2": "N/A",
        "delta_V_poly": "N/A",
        "ewald_sum": f"{ewald_multi:.2f}",
        "symmetry": "N/A"
    })
    
    # Individual Report
    report = {
        "dopant": el,
        "configurations": [
            {"type": "dilute", "file": f"{el}_oct_x0.01_dilute.xyz", "ewald": ewald_dilute, "metrics": {"avg_dist": avg_d, "lambda": lam, "sigma2": sig2}},
            {"type": "clustered", "file": f"{el}_oct_x0.01_clustered.xyz", "ewald": ewald_cluster},
            {"type": "multi", "file": f"{el}_oct_x0.05.xyz", "ewald": ewald_multi}
        ],
        "physical_reasoning": f"Octahedral site preferred for {el} if crystal field stabilization dominates (d3 high spin). VdW site minimizes steric strain (delta_V)."
    }
    with open(os.path.join(out_dir, f"{el}_oct_report.json"), "w") as f:
        json.dump(report, f, indent=4)

# Master CSV
csv_file = os.path.join(base_output_dir, "octahedral_summary.csv")
csv_columns = ["dopant", "concentration", "config_type", "avg_d_nn", "lambda", "sigma2", "delta_V_poly", "ewald_sum", "symmetry"]

with open(csv_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=csv_columns, delimiter="|")
    writer.writeheader()
    for row in summary_rows:
        writer.writerow(row)

print(f"Summary written to {csv_file}")
print("Done.")
