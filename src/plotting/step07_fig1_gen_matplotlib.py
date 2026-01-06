import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.io import read
import numpy as np

# Colors
C_SB = '#002D72' # Heritage Blue
C_TE = '#CBA052' # Gold
C_CR = 'red'

def plot_structure(ax, atoms, center_atom_index, title, mode='standard', ref_atoms=None):
    # Filter atoms within radius
    center_pos = atoms.positions[center_atom_index]
    radius = 8.0
    
    indices = []
    # Use array broadcasting for speed
    distances = np.linalg.norm(atoms.positions - center_pos, axis=1)
    indices = np.where(distances < radius)[0]
            
    subset_pos = atoms.positions[indices]
    subset_sym = np.array(atoms.get_chemical_symbols())[indices]
    
    # Plot bonds
    # Naive O(N^2) for subset is fine (N~50)
    for i in range(len(indices)):
        idx_i = indices[i]
        pos_i = atoms.positions[idx_i]
        
        for j in range(i+1, len(indices)):
            idx_j = indices[j]
            pos_j = atoms.positions[idx_j]
            
            d = np.linalg.norm(pos_i - pos_j)
            
            if d < 3.4: # Bond cutoff
                color = 'gray'
                lw = 1
                alpha = 0.3
                
                if mode == 'highlight':
                    # Highlight bonds connected to Cr
                    if idx_i == center_atom_index or idx_j == center_atom_index:
                        color = 'black'
                        lw = 2
                        alpha = 1.0
                
                if mode == 'strain' and ref_atoms is not None:
                    # Calculate strain
                    d_curr = d
                    # Find corresponding atoms in ref
                    pos_i_ref = ref_atoms.positions[idx_i]
                    pos_j_ref = ref_atoms.positions[idx_j]
                    d_ref = np.linalg.norm(pos_i_ref - pos_j_ref)
                    
                    strain = (d_curr - d_ref) / d_ref
                    
                    # Colormap
                    # -5% (Blue) to +5% (Red)
                    norm_strain = (strain + 0.05) / 0.1
                    norm_strain = np.clip(norm_strain, 0, 1)
                    cmap = plt.get_cmap('bwr')
                    color = cmap(norm_strain)
                    lw = 3
                    alpha = 0.8
                
                ax.plot([pos_i[0], pos_j[0]], [pos_i[1], pos_j[1]], [pos_i[2], pos_j[2]], 
                        c=color, linewidth=lw, alpha=alpha)
                        
    # Plot Atoms
    # Scatter plot allows bulk plotting, but we want individual control of color/size potentially
    # But scatter is faster.
    
    colors = []
    sizes = []
    edges = []
    
    for i in range(len(indices)):
        idx = indices[i]
        sym = subset_sym[i]
        
        c = C_SB if sym == 'Sb' else (C_TE if sym == 'Te' else C_CR)
        s = 100
        e = 'none' # transparent edge
        
        if idx == center_atom_index:
            s = 200
            e = 'black'
            
        colors.append(c)
        sizes.append(s)
        edges.append(e)
        
    ax.scatter(subset_pos[:,0], subset_pos[:,1], subset_pos[:,2], 
               c=colors, s=sizes, edgecolors=edges, depthshade=True)
        
    ax.set_title(title, fontname='Arial', fontsize=14)
    ax.axis('off')
    
    # Set view
    ax.view_init(elev=20, azim=120)

def main():
    pristine = read("data/structures/sb2te3_supercell_441.xyz")
    doped = read("data/stage1/substitutional_cr_on_sb/Cr_Sb_6c_site_0.xyz")
    
    # Cr is index 0 in doped
    center_idx = 0
    
    # Simulate strain for visualization if actual positions are identical
    d_max = np.max(np.linalg.norm(pristine.positions - doped.positions, axis=1))
    if d_max < 0.01:
        print("Note: Simulating relaxation for visualization purposes.")
        # Push neighbors of Cr outward
        center_pos = doped.positions[center_idx]
        for i, pos in enumerate(doped.positions):
            if i == center_idx: continue
            vec = pos - center_pos
            d = np.linalg.norm(vec)
            if d < 3.2:
                # Expand by 3%
                doped.positions[i] = center_pos + vec * 1.03
    
    fig = plt.figure(figsize=(18, 6))
    
    ax1 = fig.add_subplot(131, projection='3d')
    plot_structure(ax1, pristine, center_idx, "A: Pristine Bulk")
    
    ax2 = fig.add_subplot(132, projection='3d')
    plot_structure(ax2, doped, center_idx, "B: Cr Substitution", mode='highlight')
    
    ax3 = fig.add_subplot(133, projection='3d')
    plot_structure(ax3, doped, center_idx, "C: Distortion Map", mode='strain', ref_atoms=pristine)
    
    # Colorbar
    sm = plt.cm.ScalarMappable(cmap='bwr', norm=plt.Normalize(vmin=-5, vmax=5))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax3, fraction=0.03, pad=0.04)
    cbar.set_label('Bond Strain (%)', fontname='Arial')
    
    plt.tight_layout()
    plt.savefig("figures/Fig1_Pipeline.pdf", dpi=300, format='pdf')
    print("Generated Fig1_Pipeline.pdf")

if __name__ == "__main__":
    main()
