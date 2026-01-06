import sys
import os
import numpy as np
from ase.io import read, write
from ase.build import surface, bulk as build_bulk
from ase.geometry import get_distances

# Configuration
input_cif = 'stepA_bulk_supercells/bulk_Ti_DFT.cif'
output_dir = 'data/structures/ti_twin/step1_5811_single_twin'
output_cif = os.path.join(output_dir, 'Ti_5811_twin.cif')
output_xyz = os.path.join(output_dir, 'Ti_5811_twin.extxyz')
log_file = os.path.join(output_dir, 'logs/step1_5811.log')

def get_layers(atoms, axis=2, tolerance=0.1):
    positions = atoms.get_positions()
    coords = positions[:, axis]
    argsort = np.argsort(coords)
    sorted_coords = coords[argsort]
    unique_values = [sorted_coords[0]]
    for z in sorted_coords[1:]:
        if abs(z - unique_values[-1]) > tolerance:
            unique_values.append(z)
    return np.array(unique_values)

def main():
    if not os.path.exists(os.path.dirname(log_file)):
        os.makedirs(os.path.dirname(log_file))
        
    log_lines = []
    log_lines.append("Step 1: (5 8 11) Single Twin Supercell Construction")
    log_lines.append("="*50)
    
    # 1. Load Bulk and Extract Lattice Constant
    if not os.path.exists(input_cif):
        msg = f"Error: Input file {input_cif} not found."
        print(msg)
        return

    loaded_bulk = read(input_cif)
    vol = loaded_bulk.get_volume()
    natoms = len(loaded_bulk)
    vol_per_atom = vol / natoms
    a0 = (2 * vol_per_atom)**(1/3)
    
    log_lines.append(f"Loaded bulk structure: {natoms} atoms")
    log_lines.append(f"Extracted lattice constant a0 = {a0:.5f} A")
    
    prim_bulk = build_bulk('Ti', 'bcc', a=a0)
    
    # 2. Reorient and Create Slab
    miller = (5, 8, 11)
    log_lines.append(f"\nTarget Plane: {miller}")
    
    # Create minimal surface unit cell from primitive
    # This ensures minimal lateral vectors.
    slab_unit = surface(prim_bulk, miller, layers=1, vacuum=0)
    
    log_lines.append(f"Surface Unit Cell (Minimal):\n{slab_unit.cell}")
    
    # Lateral Expansion: Ensure minimal periodic width > 15.0 A
    # Calculate widths (distances between parallel planes defined by other vectors)
    # Area_xy = |v0 x v1|
    # width_0 = Area / |v1|
    # width_1 = Area / |v0|
    # Note: This is for 2D.
    
    def get_lateral_widths(cell):
        # cell[0] and cell[1] are the lateral vectors.
        v0 = cell[0]
        v1 = cell[1]
        area = np.linalg.norm(np.cross(v0, v1))
        w0 = area / np.linalg.norm(v1) # width along direction 0 (perp to v1)
        w1 = area / np.linalg.norm(v0) # width along direction 1 (perp to v0)
        return w0, w1

    w0, w1 = get_lateral_widths(slab_unit.cell)
    log_lines.append(f"Initial lateral widths: w0={w0:.2f} A, w1={w1:.2f} A")
    
    target_width = 15.0
    rep_x = 1
    rep_y = 1
    
    if w0 < target_width:
        rep_x = int(np.ceil(target_width / w0))
    if w1 < target_width:
        rep_y = int(np.ceil(target_width / w1))
        
    if rep_x > 1 or rep_y > 1:
        log_lines.append(f"Expanding lateral cell by {rep_x}x{rep_y} to meet > {target_width} A width.")
        # Ensure PBC is set for repeat
        slab_unit.set_pbc([True, True, True])
        slab_unit = slab_unit.repeat((rep_x, rep_y, 1))
        log_lines.append(f"Expanded Surface Cell:\n{slab_unit.cell}")
    
    # 3. Construct exactly 20 atomic layers
    unique_z = get_layers(slab_unit, axis=2, tolerance=0.1)
    n_layers_per_unit = len(unique_z)
    log_lines.append(f"Atomic layers per surface unit cell: {n_layers_per_unit}")
    
    # Calculate how many surface 'layers' (units) we need to get >= 20 atomic layers
    if n_layers_per_unit == 0:
        # Should not happen
        print("Error: 0 layers in unit slab")
        return

    needed_surface_layers = int(np.ceil(20 / n_layers_per_unit))
    
    # Add a buffer layer just in case, though usually exact is fine
    # If 20 layers are needed and unit has 10, needed=2. 2*10=20.
    
    log_lines.append(f"Requesting {needed_surface_layers} surface units to reach 20 atomic layers.")
    
    # Regenerate slab with required thickness
    slab = surface(prim_bulk, miller, layers=needed_surface_layers, vacuum=0)
    
    # Apply the lateral expansion
    if rep_x > 1 or rep_y > 1:
        log_lines.append(f"Applying lateral repeat ({rep_x}, {rep_y}, 1) to final slab.")
        slab.set_pbc([True, True, True])
        slab = slab.repeat((rep_x, rep_y, 1))
    
    # Prune to exactly 20 layers
    # We want the bottom 20 layers.
    
    unique_z_super = get_layers(slab, axis=2, tolerance=0.1)
    
    if len(unique_z_super) < 20:
        log_lines.append(f"Error: Generated slab has only {len(unique_z_super)} layers, wanted 20.")
        print("Not enough layers")
        return
            
    max_z_keep = unique_z_super[19] + 0.05
    initial_slab = slab[slab.positions[:, 2] < max_z_keep]
    log_lines.append(f"Constructed initial slab with {len(initial_slab)} atoms (20 layers).")
    
    # 4. Identify Geometric Mid-Plane & 5. Reflect
    unique_z_20 = get_layers(initial_slab, axis=2, tolerance=0.1)
    z_layer_10 = unique_z_20[9] 
    
    log_lines.append(f"Twin Plane identified at Layer 10: Z = {z_layer_10:.4f} A")
    
    bottom_half = initial_slab[initial_slab.positions[:, 2] <= (z_layer_10 + 0.05)]
    
    top_half = bottom_half.copy()
    pos = top_half.get_positions()
    pos[:, 2] = 2 * z_layer_10 - pos[:, 2]
    top_half.set_positions(pos)
    
    # 6. Merge
    full_slab = bottom_half + top_half
    
    # 7. Remove Duplicates
    temp_cell = initial_slab.cell.copy()
    z_span = np.max(full_slab.positions[:, 2]) - np.min(full_slab.positions[:, 2])
    temp_cell[2, 2] = z_span + 20.0 
    full_slab.set_cell(temp_cell)
    full_slab.set_pbc([True, True, False])
    
    dists = get_distances(full_slab.positions, cell=temp_cell, pbc=[True, True, False])
    dist_mat = dists[1]
    np.fill_diagonal(dist_mat, 100.0) 
    
    tolerance = 0.05
    rows, cols = np.where((dist_mat < tolerance) & (np.triu(np.ones(dist_mat.shape), k=1).astype(bool)))
    
    duplicates = sorted(list(set(cols)), reverse=True) 
    
    if len(duplicates) > 0:
        del full_slab[duplicates]
        
    log_lines.append(f"Removed {len(duplicates)} duplicate atoms at interface.")
    log_lines.append(f"Final atom count: {len(full_slab)}")
    
    # 8. Final PBC and Vacuum
    vacuum = 15.0 
    full_slab.center(vacuum=vacuum/2.0, axis=2)
    full_slab.set_pbc([True, True, True])
    
    log_lines.append(f"Final Cell with Vacuum:\n{full_slab.cell}")
    
    write(output_cif, full_slab)
    write(output_xyz, full_slab)
    log_lines.append(f"Saved {output_cif}")
    log_lines.append(f"Saved {output_xyz}")
    
    with open(log_file, 'w') as f:
        f.write('\n'.join(log_lines))
    print("Processing complete.")

if __name__ == "__main__":
    main()