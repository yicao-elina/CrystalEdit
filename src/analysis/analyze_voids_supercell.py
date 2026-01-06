import numpy as np
from pymatgen.core import Structure
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.analysis.local_env import CrystalNN
import json
import os

from ase.io import read as ase_read
from pymatgen.io.ase import AseAtomsAdaptor

# 1. Load Supercell Structure
filename = "data/structures/sb2te3_supercell_441.xyz"
if not os.path.exists(filename):
    filename = "sb2te3_supercell_441.xyz"

print(f"Loading supercell from {filename} using ASE...")
atoms = ase_read(filename)
structure = AseAtomsAdaptor.get_structure(atoms)


# 2. Identify Interstitial Sites (Voids) using Voronoi on the Supercell
# This will find all voids in the 4x4x1 supercell
print("Identifying voids using Voronoi tessellation on supercell...")
generator = VoronoiInterstitialGenerator()
# We use a dummy species "Cr" to represent the center of the void
interstitials = list(generator.generate(structure, {"Cr": 1}))

print(f"Found {len(interstitials)} void candidates in the supercell.")

# 3. Analyze each site
cnn = CrystalNN()
voids_results = []

# To identify VdW gaps, we look for sites with large distances to neighbors 
# and specific coordination.
# In Sb2Te3, octahedral sites exist within the QL, and voids exist in the VdW gap.

for i, defect in enumerate(interstitials):
    site = defect.site
    frac_coords = site.frac_coords.tolist()
    
    # Analyze environment in the supercell
    temp_struct = structure.copy()
    temp_struct.append("Cr", frac_coords)
    
    try:
        cn_info = cnn.get_nn_info(temp_struct, len(temp_struct)-1)
        coordination_number = len(cn_info)
        
        neighbor_data = []
        for neighbor in cn_info:
            n_site = neighbor['site']
            dist = site.distance(n_site)
            neighbor_data.append({
                "element": n_site.specie.symbol,
                "distance": float(dist)
            })
            
        # Sort neighbors by distance
        neighbor_data = sorted(neighbor_data, key=lambda x: x['distance'])
        min_dist = neighbor_data[0]['distance'] if neighbor_data else 0.0
        
        # Classification
        # Octahedral: CN=6
        # Tetrahedral: CN=4
        # VdW gap: Usually sites between Te layers with larger distances 
        # or specific Z-planes.
        
        # Logic for classification:
        if coordination_number == 4:
            void_type = "tet"
        elif coordination_number == 6:
            void_type = "oct"
        elif min_dist > 2.8: # Heuristic for VdW gap sites being "looser"
            void_type = "vdw"
        else:
            # Check coordination geometry for distorted oct/tet
            if coordination_number > 6:
                void_type = "vdw" # Often VdW sites have high effective coordination
            else:
                void_type = f"distorted_{coordination_number}"

        voids_results.append({
            "void_type": void_type,
            "frac_coords": frac_coords,
            "nearest_neighbors": neighbor_data
        })
        
    except Exception as e:
        print(f"Error analyzing void {i}: {e}")

# 4. Filter or label the VdW gaps more precisely if possible
# The largest interlayer spacing in Sb2Te3 is the Te-Te gap.
# In conventional cell z-coords, these are at approx 0.11, 0.44, 0.78 (centers)
# In supercell (4x4x1), z-coords are same.

# Sort results by Z for better organization
voids_results = sorted(voids_results, key=lambda x: x['frac_coords'][2])

# 5. Output JSON
output_json = "void_space_analysis_supercell.json"
if os.path.exists("26ICML-CrystalEdit"):
    output_json = "data/reports/void_space_analysis_supercell.json"

with open(output_json, "w") as f:
    json.dump(voids_results, f, indent=4)

print(f"Void space analysis completed. Results saved to {output_json}")
print(f"Summary of supercell voids:")
types = [v['void_type'] for v in voids_results]
from collections import Counter
counts = collections.Counter(types) if 'collections' in globals() else Counter(types)
for t, c in counts.items():
    print(f"  {t}: {c}")
