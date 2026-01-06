import numpy as np
from pymatgen.core import Structure
from pymatgen.analysis.defects.generators import VoronoiInterstitialGenerator
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import os

# 1. Load Structure
filename = "Sb2Te3-mp-1201.cif"
if not os.path.exists(filename):
    print(f"File {filename} not found in CWD. Trying with dir prefix...")
    filename = "data/structures/Sb2Te3-mp-1201.cif"

print(f"Loading structure from {filename}...")
structure = Structure.from_file(filename)

# 2. Identify Interstitial Sites
print("Generating interstitial sites...")
# We use VoronoiInterstitialGenerator. 
# clustering_tol can be adjusted if we want to distinguish sites that are very close.
# Default is usually around 0.5 Angstrom. Let's stick to default to find physically distinct ones.
generator = VoronoiInterstitialGenerator()
interstitials = list(generator.generate(structure, {"Cr": 1}))

print(f"Found {len(interstitials)} distinct interstitial sites candidates.")

# Prepare outputs
site_reasoning_md = "# Interstitial Site Analysis & Reasoning\n\n"
site_reasoning_md += "Based on Voronoi tessellation and local environment analysis.\n\n"
site_reasoning_md += "## Methodology\n"
site_reasoning_md += "1.  **Generation**: Voronoi tessellation was used to identify the largest voids in the lattice compatible with Cr insertion.\n"
site_reasoning_md += "2.  **Analysis**: For each distinct site, `CrystalNN` was used to determine the coordination environment and nearest neighbors.\n"
site_reasoning_md += "3.  **Classification**: Sites are classified by their coordination number and the chemical identity of the neighbors.\n\n"
site_reasoning_md += "## Identified Sites\n"

neighbor_distances_csv = [["Site_Index", "Label", "Coordination_Number", "Nearest_Neighbor", "Distance_A", "Environment_Description"]]

xyz_coordinates_for_plot = []
xyz_labels_for_plot = []

# Supercell scaling matrix
scaling_matrix = [[4, 0, 0], [0, 4, 0], [0, 0, 1]]

# Local Environment Analyzer
cnn = CrystalNN()

for i, defect in enumerate(interstitials):
    site = defect.site
    frac_coords = site.frac_coords
    cart_coords = site.coords
    
    # Create a temporary structure with the interstitial to analyze local env
    temp_structure = structure.copy()
    temp_structure.append("Cr", frac_coords)
    
    # Analyze the inserted atom (last index)
    try:
        cn_info = cnn.get_nn_info(temp_structure, len(temp_structure)-1)
        coordination_number = len(cn_info)
        
        # Analyze neighbors
        neighbor_elements = []
        distances = []
        for neighbor in cn_info:
            site_obj = neighbor['site']
            neighbor_elements.append(site_obj.specie.symbol)
            # Distance is implicitly handled by CrystalNN weights, but let's get actual distance
            distances.append(site.distance(site_obj))
            
        min_dist = min(distances) if distances else 0.0
        nearest_neighbor_el = neighbor_elements[distances.index(min_dist)] if distances else "None"
        
        # Categorize
        neighbor_counts = {el: neighbor_elements.count(el) for el in set(neighbor_elements)}
        env_desc = f"CN{coordination_number}-" + "-".join([f"{k}{v}" for k,v in neighbor_counts.items()])
        
        # Determine likely hole type based on CN and geometry (simplified)
        if coordination_number == 6:
            polyhedra = "Octahedral-like"
        elif coordination_number == 4:
            polyhedra = "Tetrahedral-like"
        else:
            polyhedra = f"Distorted-{coordination_number}-coord"
            
    except Exception as e:
        print(f"Error analyzing site {i}: {e}")
        coordination_number = "Unknown"
        env_desc = "Unknown"
        min_dist = 0.0
        nearest_neighbor_el = "None"
        polyhedra = "Unknown"

    # Label for file and plot
    label = f"{polyhedra}_{env_desc}"
    clean_label = label.replace(" ", "_")
    
    # 1. Create 4x4x1 Supercell with defect
    sc_structure = structure.copy()
    sc_structure.make_supercell(scaling_matrix)
    
    # Map coordinates to supercell
    sc_frac_coords = [frac_coords[0]/4.0, frac_coords[1]/4.0, frac_coords[2]/1.0]
    sc_structure.append("Cr", sc_frac_coords)
    
    # Reason String
    reason = (
              f"Voronoi-generated interstitial site #{i}. "
              f"Geometry: {polyhedra}. "
              f"Coordination: {coordination_number} ({env_desc}). "
              f"Nearest Neighbor: {nearest_neighbor_el} at {min_dist:.3f} A."
             )
    
    # Append to data
    neighbor_distances_csv.append([i, clean_label, coordination_number, nearest_neighbor_el, f"{min_dist:.4f}", env_desc])
    xyz_coordinates_for_plot.append(cart_coords)
    xyz_labels_for_plot.append(f"{i}:{env_desc}")
    
    # Update MD
    site_reasoning_md += f"### Site {i}: {polyhedra}\n"
    site_reasoning_md += f"- **Label**: `{clean_label}`\n"
    site_reasoning_md += f"- **Coordinates (Frac)**: ({frac_coords[0]:.3f}, {frac_coords[1]:.3f}, {frac_coords[2]:.3f})\n"
    site_reasoning_md += f"- **Coordination**: {coordination_number}\n"
    site_reasoning_md += f"- **Environment**: {env_desc} (Neighbors: {neighbor_elements})\n"
    site_reasoning_md += f"- **Nearest Neighbor Distance**: {min_dist:.3f} Å\n"
    site_reasoning_md += f"- **Reason for Selection**: Identified as a distinct Voronoi node representing a potential local minimum for intercalation. {polyhedra} geometry suggests stability for transition metals.\n\n"
    
    # Write Extended XYZ
    atoms = AseAtomsAdaptor.get_atoms(sc_structure)
    # Add reason to atoms info? ASE extxyz writer handles info dict
    atoms.info['description'] = reason
    atoms.info['site_type'] = clean_label
    
    filename = f"Cr_interstitial_{i}_{clean_label}.xyz"
    # Ensure filename is safe
    filename = filename.replace("/", "_")
    
    # If running in root, prefix with folder
    if not os.path.exists("26ICML-CrystalEdit") and os.path.exists("../26ICML-CrystalEdit"):
         filename = f"../26ICML-CrystalEdit/{filename}"
    elif os.path.exists("26ICML-CrystalEdit"):
         filename = f"26ICML-CrystalEdit/{filename}"
         
    write(filename, atoms, format="extxyz")
    print(f"Written {filename}")

# Save CSV
csv_filename = "neighbor_distances.csv"
if os.path.exists("26ICML-CrystalEdit"):
    csv_filename = "data/tables/neighbor_distances.csv"
    
with open(csv_filename, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(neighbor_distances_csv)

# Save Markdown
md_filename = "site_reasoning.md"
if os.path.exists("26ICML-CrystalEdit"):
    md_filename = "docs/site_reasoning.md"

with open(md_filename, "w") as f:
    f.write(site_reasoning_md)

# Visualization Plot
print("Generating visualization...")
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Host atoms
host_coords = structure.cart_coords
host_species = [s.symbol for s in structure.species]
colors = {'Sb': 'blue', 'Te': 'orange'}

# Scatter host
for j, coord in enumerate(host_coords):
    el = host_species[j]
    ax.scatter(coord[0], coord[1], coord[2], c=colors[el], alpha=0.15, s=30) # lighter alpha to see inside

# Dummy points for legend
ax.scatter([], [], [], c='blue', label='Sb', alpha=1)
ax.scatter([], [], [], c='orange', label='Te', alpha=1)

# Plot interstitials
xyz_coordinates_for_plot = np.array(xyz_coordinates_for_plot)
if len(xyz_coordinates_for_plot) > 0:
    # Use a colormap for different sites
    cmap = plt.get_cmap('Set1')
    
    for i, (coord, label) in enumerate(zip(xyz_coordinates_for_plot, xyz_labels_for_plot)):
        ax.scatter(coord[0], coord[1], coord[2], color=cmap(i%9), s=150, marker='*', label=f"Site {i}")
        ax.text(coord[0], coord[1], coord[2], f" {i}", color='black', fontsize=12, fontweight='bold')

ax.set_xlabel('X (Å)')
ax.set_ylabel('Y (Å)')
ax.set_zlabel('Z (Å)')
ax.set_title(r'$Sb_2Te_3$ Interstitial Sites Locations')
ax.legend(loc='upper right')

pdf_filename = "interstitial_sites_3d.pdf"
if os.path.exists("26ICML-CrystalEdit"):
    pdf_filename = "figures/interstitial_sites_3d.pdf"

plt.savefig(pdf_filename)
print(f"Written visualization to {pdf_filename}.")
print("Done.")
