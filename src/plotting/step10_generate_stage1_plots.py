import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.io import read
import numpy as np
import os
from scipy.spatial import cKDTree

# JHU Colors
C_SB = '#002D72' # Heritage Blue
C_TE = '#CBA052' # Gold
C_DEFECT = '#FF0000' # Spirit Red / Highlight

def get_color(symbol):
    if symbol == 'Sb': return C_SB
    if symbol == 'Te': return C_TE
    return C_DEFECT # H, Cr, V, Mn

def create_legend(ax, symbols):
    from matplotlib.lines import Line2D
    elements = sorted(list(set(symbols)))
    legend_elements = []
    for el in elements:
        c = get_color(el)
        legend_elements.append(Line2D([0], [0], marker='o', color='w', label=el, markerfacecolor=c, markersize=10))
    
    legend_elements.append(Line2D([0], [0], color='blue', lw=2, label='Compressive'))
    legend_elements.append(Line2D([0], [0], color='red', lw=2, label='Tensile'))
    
    ax.legend(handles=legend_elements, loc='upper right', title="Species")

def plot_structure(ax, atoms, title, mode='standard', ref_atoms=None):
    pos = atoms.positions
    sym = atoms.get_chemical_symbols()
    
    # Crop for visualization? 
    # Supercell is large. Let's focus on the center of the cell or defect.
    # Identify "Defect Center": 
    # - If interstitial: the new atom (last usually)
    # - If substitution: the atom with changed symbol
    # - If vacancy: hard to define from XYZ alone without metadata.
    # Strategy: Plot everything but use depth shading.
    
    # 1. Strain Mapping Setup
    mapping = {} # def_idx -> ref_idx
    if mode == 'strain' and ref_atoms is not None:
        tree = cKDTree(ref_atoms.positions)
        dists, indices = tree.query(pos, k=1)
        for i, (d, ref_idx) in enumerate(zip(dists, indices)):
            if d < 0.8: # Threshold for mapping (relaxed atoms move < 0.8 A)
                mapping[i] = ref_idx

    # 2. Plot Bonds
    # Using cKDTree for neighbor finding
    tree_self = cKDTree(pos)
    pairs = tree_self.query_pairs(r=3.4)
    
    for i, j in pairs:
        p1 = pos[i]
        p2 = pos[j]
        
        col = 'gray'
        lw = 0.5
        alpha = 0.3
        
        if mode == 'strain' and ref_atoms is not None:
            if i in mapping and j in mapping:
                ref_i = mapping[i]
                ref_j = mapping[j]
                
                # Check if bond existed in reference
                p1_ref = ref_atoms.positions[ref_i]
                p2_ref = ref_atoms.positions[ref_j]
                d_ref = np.linalg.norm(p1_ref - p2_ref)
                
                if d_ref < 3.5: # Valid ref bond
                    d_curr = np.linalg.norm(p1 - p2)
                    strain = (d_curr - d_ref) / d_ref
                    
                    # Map strain
                    norm_strain = (strain + 0.05) / 0.1
                    norm_strain = np.clip(norm_strain, 0, 1)
                    col = plt.get_cmap('bwr')(norm_strain)
                    lw = 1.5
                    alpha = 0.8
        
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], c=col, lw=lw, alpha=alpha)

    # 3. Plot Atoms
    colors = [get_color(s) for s in sym]
    sizes = [60] * len(atoms)
    
    # Highlight non-host
    for i, s in enumerate(sym):
        if s not in ['Sb', 'Te']:
            sizes[i] = 150 # Larger for dopants
            
    ax.scatter(pos[:,0], pos[:,1], pos[:,2], c=colors, s=sizes, depthshade=True, edgecolors='black', linewidth=0.5)
    
    ax.set_title(title, fontname='Arial', fontsize=12)
    ax.axis('off')
    ax.view_init(elev=20, azim=120)

def generate_pdf(xyz_path, ref_path, output_path):
    try:
        atoms = read(xyz_path)
        ref = read(ref_path)
    except Exception as e:
        print(f"Error reading {xyz_path}: {e}")
        return

    name = os.path.basename(xyz_path).replace('.xyz', '')
    
    fig = plt.figure(figsize=(18, 6))
    
    # Panel A: Reference
    ax1 = fig.add_subplot(131, projection='3d')
    plot_structure(ax1, ref, "A: Pristine Host")
    
    # Panel B: Defect Structure
    ax2 = fig.add_subplot(132, projection='3d')
    plot_structure(ax2, atoms, f"B: {name}")
    
    # Panel C: Strain Map
    ax3 = fig.add_subplot(133, projection='3d')
    plot_structure(ax3, atoms, "C: Bond Strain Map", mode='strain', ref_atoms=ref)
    
    # Legend
    create_legend(ax3, atoms.get_chemical_symbols())
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, format='pdf')
    plt.close()
    print(f"Generated {output_path}")

def main():
    root_dir = "data/stage1"
    ref_file = "data/structures/sb2te3_supercell_441.xyz"
    
    if not os.path.exists(ref_file):
        print("Reference file not found!")
        return

    # Walk through stage1
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for f in filenames:
            if f.endswith(".xyz"):
                full_path = os.path.join(dirpath, f)
                
                # Determine output name
                # Structure: stage1/category/file.xyz -> stage1/category/Fig_file.pdf
                out_name = f"Fig_{f.replace('.xyz', '.pdf')}"
                out_path = os.path.join(dirpath, out_name)
                
                generate_pdf(full_path, ref_file, out_path)

if __name__ == "__main__":
    main()
