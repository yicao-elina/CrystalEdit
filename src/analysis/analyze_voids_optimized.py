import numpy as np
from pymatgen.core import Structure
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.analysis.local_env import CrystalNN
from ase.io import read as ase_read
from pymatgen.io.ase import AseAtomsAdaptor
import json
import os

# 1. Load Host Unit Cell
cif_filename = "data/structures/Sb2Te3-mp-1201.cif"
if not os.path.exists(cif_filename):
    cif_filename = "Sb2Te3-mp-1201.cif"

print(f"Loading unit cell from {cif_filename}...")
unit_cell = Structure.from_file(cif_filename)

# 2. Identify VdW Gaps
print("Identifying Van der Waals gap locations...")
# Sort z-coordinates of planes
z_coords = np.unique(np.round(unit_cell.frac_coords[:, 2], 4))
z_coords = np.sort(z_coords)
# Include periodic boundary
z_extended = np.append(z_coords, z_coords[0] + 1.0)
diffs = np.diff(z_extended)
max_gap_idx = np.argmax(diffs)
vdw_gap_center_z = (z_extended[max_gap_idx] + z_extended[max_gap_idx+1]) / 2.0
if vdw_gap_center_z > 1.0: vdw_gap_center_z -= 1.0

print(f"VdW gap center identified at z = {vdw_gap_center_z:.4f}")

# 3. Identify and Classify Voids in Unit Cell
print("Identifying unique voids in unit cell...")
generator = VoronoiInterstitialGenerator()
unit_interstitials = list(generator.generate(unit_cell, {"Cr": 1}))

cnn = CrystalNN()
classified_unit_voids = []

for defect in unit_interstitials:
    site = defect.site
    f_coords = site.frac_coords
    
    temp_unit = unit_cell.copy()
    temp_unit.append("Cr", f_coords)
    
    try:
        cn_info = cnn.get_nn_info(temp_unit, len(temp_unit)-1)
        cn = len(cn_info)
        
        neighbor_data = []
        for neighbor in cn_info:
            n_site = neighbor['site']
            dist = site.distance(n_site)
            neighbor_data.append({"element": n_site.specie.symbol, "distance": float(dist)})
        
        neighbor_data = sorted(neighbor_data, key=lambda x: x['distance'])
        distances = [x['distance'] for x in neighbor_data]
        
        # Check VdW Proximity
        z_diff = abs(f_coords[2] - vdw_gap_center_z)
        z_diff = min(z_diff, 1.0 - z_diff)
        is_vdw_loc = z_diff < 0.08
        
        # Initialize flags
        is_tet = False
        is_oct = False
        core_cn = cn
        
        # Determine Geometry Type
        # Check for Tetrahedral (4 close)
        if len(distances) >= 4:
            if len(distances) == 4 or (distances[4] - distances[3] > 0.2):
                is_tet = True
                core_cn = 4
        
        # Check for Octahedral (6 close)
        if len(distances) >= 6:
            # 6 neighbors within a reasonable shell
            if len(distances) == 6 or (distances[6] - distances[5] > 0.2):
                is_oct = True
                core_cn = 6
            elif cn == 7: # If CrystalNN says 7, it's likely a distorted 6+1
                is_oct = True
                core_cn = 6
            elif cn == 8 and is_vdw_loc: 
                is_oct = True 
                core_cn = 6

        min_dist = distances[0] if distances else 0.0
        
        # Final Naming
        if is_tet:
            void_type = "Tetrahedral"
        elif is_oct:
            void_type = "Octahedral"
        elif cn == 8:
            # Explicitly map CN=8 to VdW Octahedral for this structure
            void_type = "Octahedral"
            is_vdw_loc = True # Force VdW label
        else:
            void_type = f"Distorted (CN={cn})"
            
        # Add Location Context
        if is_vdw_loc:
            void_type += " (VdW Gap)"
        else:
            void_type += " (Intralayer)"
            
        classified_unit_voids.append({
            "void_type": void_type,
            "frac_coords_unit": f_coords.tolist(),
            "nn_unit": neighbor_data
        })
    except Exception as e:
        print(f"Error classifying: {e}")
        pass

# 4. Expand to Supercell (4x4x1)
print("Mapping voids to 4x4x1 supercell...")
supercell_filename = "data/structures/sb2te3_supercell_441.xyz"
if not os.path.exists(supercell_filename):
    supercell_filename = "sb2te3_supercell_441.xyz"

atoms_sc = ase_read(supercell_filename)
structure_sc = AseAtomsAdaptor.get_structure(atoms_sc)

final_voids = []

for uv in classified_unit_voids:
    f_unit = uv['frac_coords_unit']
    for i in range(4):
        for j in range(4):
            # Supercell fractional coordinates
            # f_sc_x = (f_unit_x + i) / 4.0
            # f_sc_y = (f_unit_y + j) / 4.0
            # f_sc_z = f_unit_z / 1.0
            f_sc = [(f_unit[0] + i)/4.0, (f_unit[1] + j)/4.0, f_unit[2]]
            
            # Recalculate NN in supercell for accuracy (or just map distances)
            # Actually, distances are the same, but the JSON should be for the supercell.
            
            # We need elements and distances for the supercell structure
            # temp_sc = structure_sc.copy()
            # temp_sc.append("Cr", f_sc)
            # site_sc = temp_sc[-1]
            # neighbors = structure_sc.get_neighbors(site_sc, r=5.0)
            # ... this might be slow for all sites. 
            # Since it's a supercell expansion, elements and distances are inherited.
            
            final_voids.append({
                "void_type": uv['void_type'],
                "frac_coords": f_sc,
                "nearest_neighbors": uv['nn_unit']
            })

# 5. Output JSON
output_path = "void_analysis_supercell.json"
if os.path.exists("26ICML-CrystalEdit"):
    output_path = "data/reports/void_analysis_supercell.json"

with open(output_path, "w") as f:
    json.dump(final_voids, f, indent=4)

print(f"Analysis saved to {output_path}")
print(f"Total sites in supercell: {len(final_voids)}")
