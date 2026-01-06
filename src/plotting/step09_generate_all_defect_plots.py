import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.io import read
import numpy as np
import os

# JHU Colors
C_SB = '#002D72' # Heritage Blue
C_TE = '#CBA052' # Gold
C_DEFECT = 'red' # Highlight color

def get_color(symbol):
    if symbol == 'Sb': return C_SB
    if symbol == 'Te': return C_TE
    return 'gray'

def create_legend(ax):
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Sb', markerfacecolor=C_SB, markersize=10),
        Line2D([0], [0], marker='o', color='w', label='Te', markerfacecolor=C_TE, markersize=10),
        Line2D([0], [0], color='blue', lw=2, label='Compressive'),
        Line2D([0], [0], color='red', lw=2, label='Tensile')
    ]
    ax.legend(handles=legend_elements, loc='upper right', title="Legend")

def plot_atoms_and_bonds(ax, atoms, title, mode='standard', view='perspective', ref_atoms=None, highlight_z=None):
    # Filter for visualization speed/clarity?
    # For large supercells (GB ~400 atoms), plotting all is heavy but feasible with scatter.
    # We might want to crop to a region of interest for GB/Twin.
    
    pos = atoms.positions
    sym = atoms.get_chemical_symbols()
    
    # 1. Determine subsets based on view/highlight
    mask = np.ones(len(atoms), dtype=bool)
    
    # Crop for clarity if structure is very large? 
    # Let's keep full structure but adjust camera.
    
    colors = [get_color(s) for s in sym]
    sizes = [60] * len(atoms)
    
    # Highlight logic (e.g. interface atoms)
    if highlight_z is not None:
        # Find atoms near Z plane
        z_coords = pos[:, 2]
        # Wrap z if pbc? Assumed cartesian.
        for i, z in enumerate(z_coords):
            if abs(z - highlight_z) < 2.5:
                # Highlight edges?
                pass

    # Bonds
    # Calculate pairwise distances (naive N^2 is OK for N<500)
    # Using simple loop to allow strain coloring
    
    # Optimize: Pre-calculate neighbor list?
    # Use ASE neighborlist if available, else brute force
    from ase.neighborlist import neighbor_list
    # cutoffs: Sb-Te ~ 3.2. 
    cutoffs = [1.8] * len(atoms) # radii
    # nl = neighbor_list('ij', atoms, cutoffs) -> too standard.
    
    # Brute force bond drawing for custom coloring
    # Only draw bonds for a subset to avoid clutter?
    
    # Draw bonds
    # Strain calculation needs reference.
    # Since geometry changed (atoms removed/added/shifted), 1-to-1 mapping is hard for GB/Slab.
    # We will simulate strain visualization for "Mode=Strain" based on bond length deviation from ideal (3.0 A).
    ideal_bond = 3.0
    
    # Limit bond drawing to "front" atoms or slice?
    # Let's draw all bonds but thin.
    
    # Distance matrix
    # Only compute upper triangle
    N = len(atoms)
    
    # To speed up, use cKDTree if scipy available (it is).
    from scipy.spatial import cKDTree
    tree = cKDTree(pos)
    pairs = tree.query_pairs(r=3.4) # Cutoff 3.4
    
    for i, j in pairs:
        p1 = pos[i]
        p2 = pos[j]
        d = np.linalg.norm(p1 - p2)
        
        col = 'gray'
        lw = 0.5
        alpha = 0.3
        
        if mode == 'strain':
            strain = (d - ideal_bond) / ideal_bond
            # Map -5% to +5%
            norm_strain = (strain + 0.05) / 0.1
            norm_strain = np.clip(norm_strain, 0, 1)
            col = plt.get_cmap('bwr')(norm_strain)
            lw = 1.5
            alpha = 0.8
            
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], c=col, lw=lw, alpha=alpha)

    # Draw Atoms
    ax.scatter(pos[:,0], pos[:,1], pos[:,2], c=colors, s=sizes, depthshade=True, edgecolors='black', linewidth=0.5)
    
    ax.set_title(title, fontname='Arial', fontsize=12)
    ax.axis('off')
    
    # View angles
    if view == 'side':
        ax.view_init(elev=0, azim=0) # X-Z plane?
        # Check lattice to align.
        # Usually c is Z. 
        # Side view: Look along Y (azim=0) or X (azim=90).
        ax.view_init(elev=5, azim=10) # Slight perspective
    elif view == 'top':
        ax.view_init(elev=90, azim=-90) # X-Y plane
    else:
        ax.view_init(elev=20, azim=120)

def generate_plot(defect_name, file_path, ref_path, output_name):
    print(f"Processing {defect_name}...")
    
    try:
        atoms = read(file_path)
    except Exception as e:
        print(f"Failed to read {file_path}: {e}")
        return

    try:
        ref = read(ref_path)
    except:
        ref = atoms # Fallback

    fig = plt.figure(figsize=(18, 6))
    
    # Panel A: Reference (Bulk)
    # Or Reference Structure (Pristine)
    ax1 = fig.add_subplot(131, projection='3d')
    plot_atoms_and_bonds(ax1, ref, f"A: Pristine Bulk (Ref)", view='perspective')
    
    # Panel B: Defect View
    ax2 = fig.add_subplot(132, projection='3d')
    
    # Determine best view based on defect type
    view_mode = 'perspective'
    if 'slab' in defect_name.lower(): view_mode = 'side'
    if 'fault' in defect_name.lower(): view_mode = 'side'
    if 'twin' in defect_name.lower(): view_mode = 'side'
    if 'gb' in defect_name.lower(): view_mode = 'top' # Top view for kites
    
    plot_atoms_and_bonds(ax2, atoms, f"B: {defect_name}", view=view_mode)
    
    # Panel C: Strain / Physics
    ax3 = fig.add_subplot(133, projection='3d')
    plot_atoms_and_bonds(ax3, atoms, f"C: Strain/Bond Map", mode='strain', view=view_mode)
    
    # Add Legend to Panel C or main figure
    create_legend(ax3)
    
    plt.tight_layout()
    plt.savefig(output_name, dpi=300, format='pdf')
    plt.close()
    print(f"Saved {output_name}")

def main():
    ref_file = "data/structures/sb2te3_supercell_441.xyz"
    
    tasks = [
        ("Te-Terminated Slab", "data/structures/sb2te3_slab_5QL_Te_term.xyz", "figures/Fig_Slab_Te.pdf"),
        ("Sb-Terminated Slab", "data/structures/sb2te3_slab_5QL_Sb_term.xyz", "figures/Fig_Slab_Sb.pdf"),
        ("Stacking Fault", "data/structures/sb2te3_faulted_stacking.xyz", "figures/Fig_Fault.pdf"),
        ("Twin Boundary", "data/structures/sb2te3_twin_boundary.xyz", "figures/Fig_Twin.pdf"),
        ("Sigma7 Grain Boundary", "data/structures/sb2te3_gb.xyz", "figures/Fig_GB.pdf")
    ]
    
    for name, path, out in tasks:
        if os.path.exists(path):
            generate_plot(name, path, ref_file, out)
        else:
            print(f"Skipping {name}, file not found: {path}")

if __name__ == "__main__":
    main()
