import numpy as np
from pymatgen.core import Structure
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.analysis.local_env import CrystalNN
import json
import os

# 1. Load Structure
filename = "Sb2Te3-mp-1201.cif"
if not os.path.exists(filename):
    print(f"File {filename} not found in CWD. Trying with dir prefix...")
    filename = "data/structures/Sb2Te3-mp-1201.cif"

print(f"Loading structure from {filename}...")
structure = Structure.from_file(filename)

# 2. Identify Interstitial Sites
print("Generating interstitial sites using Voronoi tessellation...")
generator = VoronoiInterstitialGenerator()
interstitials = list(generator.generate(structure, {"Cr": 1}))

print(f"Found {len(interstitials)} distinct interstitial sites.")

# 3. Analyze each site to classify as Tetrahedral, Octahedral, or Van der Waals
cnn = CrystalNN()
voids_data = []

for i, defect in enumerate(interstitials):
    site = defect.site
    frac_coords = site.frac_coords.tolist()
    
    # Create temp structure to analyze local environment
    temp_structure = structure.copy()
    temp_structure.append("Cr", frac_coords)
    
    # Analyze coordination
    try:
        cn_info = cnn.get_nn_info(temp_structure, len(temp_structure)-1)
        coordination_number = len(cn_info)
        
        # Determine nearest neighbors
        neighbor_distances = []
        neighbor_elements = []
        for neighbor in cn_info:
            site_obj = neighbor['site']
            dist = site.distance(site_obj)
            neighbor_distances.append(dist)
            neighbor_elements.append(site_obj.specie.symbol)
            
        min_dist = min(neighbor_distances) if neighbor_distances else 0.0
        
        # Heuristic Classification
        # Octahedral: CN=6
        # Tetrahedral: CN=4
        # Van der Waals: Large distance, usually located in gap z-coordinates
        
        # Analyze fractional Z to check if in VdW gap
        # Sb2Te3 is layered. VdW gaps are around specific Z planes in the conventional cell.
        # Structure is R-3m (166). Te1-Te1 gap.
        
        # Let's rely on coordination and geometry.
        if coordination_number == 6:
            void_type = "octahedral"
        elif coordination_number == 4:
            void_type = "tetrahedral"
        else:
            # Check if it's in a large void (VdW)
            # VdW gaps often allow higher coordination or unusual geometry with larger distances
            if min_dist > 2.8: # Threshold for "large" void
                void_type = "vdw"
            else:
                void_type = f"distorted_{coordination_number}"

        nearest_neighbors_data = []
        for el, dist in zip(neighbor_elements, neighbor_distances):
            nearest_neighbors_data.append({"element": el, "distance": dist})
            
        void_entry = {
            "void_type": void_type,
            "frac_coords": frac_coords,
            "nearest_neighbors": sorted(nearest_neighbors_data, key=lambda x: x['distance'])
        }
        voids_data.append(void_entry)
        
    except Exception as e:
        print(f"Error analyzing site {i}: {e}")

# 4. Filter for specific types if duplicates exist (Voronoi returns unique sites already)
# Just organizing output.

output_data = {
    "tetrahedral_voids": [v for v in voids_data if "tetrahedral" in v['void_type']],
    "octahedral_voids": [v for v in voids_data if "octahedral" in v['void_type']],
    "vdw_gaps": [v for v in voids_data if "vdw" in v['void_type']],
    "other_voids": [v for v in voids_data if v['void_type'] not in ["tetrahedral", "octahedral", "vdw"]]
}

# 5. Write JSON
output_filename = "void_analysis.json"
if os.path.exists("26ICML-CrystalEdit"):
    output_filename = "data/reports/void_analysis.json"

with open(output_filename, "w") as f:
    json.dump(output_data, f, indent=4)

print(f"Written void analysis to {output_filename}")
print("Summary:")
print(f"Tetrahedral: {len(output_data['tetrahedral_voids'])}")
print(f"Octahedral: {len(output_data['octahedral_voids'])}")
print(f"VdW Gaps: {len(output_data['vdw_gaps'])}")
print(f"Other: {len(output_data['other_voids'])}")
