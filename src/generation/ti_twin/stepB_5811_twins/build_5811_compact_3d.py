import os
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.build import make_supercell
from ase.geometry import get_distances

# Configuration
input_cif = 'data/structures/ti_twin/Ti-bcc.cif'
output_dir = 'data/structures/ti_twin/stepB_5811_twins'
log_file = 'docs/logs/stepB_5811_twins_compact.log'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def build_compact_twin():
    log_lines = []
    log_lines.append("Building Compact 3D (5 8 11) Twin Supercell (Corrected)")
    log_lines.append("="*40)
    
    # 1. Load Bulk
    try:
        bulk = read(input_cif)
    except Exception as e:
        print(f"Error reading input: {e}")
        return

    # 2. Orientation Matrix P (rows are new basis vectors in terms of old)
    # New Basis:
    # u = [1, -2, 1]
    # v = [5, 1, -3]
    # w = [5, 8, 11] (Normal)
    P = np.array([
        [1, -2, 1],
        [5, 1, -3],
        [5, 8, 11]
    ])
    
    # Generate Oriented Cell
    # This creates the periodic unit cell defined by P
    oriented_unit = make_supercell(bulk, P)
    
    # Expand laterally (2x in X) to meet size requirements (~16 A)
    # Original X is ~8 A.
    oriented_unit *= (2, 1, 1)
    
    log_lines.append(f"Oriented Unit Cell: {oriented_unit.get_global_number_of_atoms()} atoms")
    log_lines.append(f"Cell Dimensions: {oriented_unit.cell.cellpar()[:3]}")

    # 3. Prepare for Twin Construction
    # We will construct the twin by stacking: Grain A + Grain B.
    # Grain A: Bottom half of the unit cell.
    # Grain B: Mirror of Grain A.
    
    # First, ensure all atoms are wrapped strictly into the cell
    oriented_unit.wrap()
    
    # Calculate interplanar spacing for (5 8 11)
    # d = a / sqrt(h^2+k^2+l^2)
    # But we can also get it from the geometry of the generated cell.
    # The cell Z vector corresponds to [5 8 11].
    # Number of planes in this period?
    # BCC lattice points satisfy h+k+l = even (for conventional).
    # If we use the Primitive cell, the density is different.
    # Let's infer d_spacing from the Z-coordinates of atoms.
    
    positions = oriented_unit.get_positions()
    z_coords = positions[:, 2]
    # Sort and find unique layers
    z_unique = np.unique(np.round(z_coords, 4))
    z_diffs = np.diff(z_unique)
    if len(z_diffs) > 0:
        d_spacing_measured = np.mean(z_diffs) # Approximate
        log_lines.append(f"Measured avg Z-spacing: {d_spacing_measured:.4f} A")
    else:
        d_spacing_measured = 0.2 # Fallback
        
    # 4. Select Grain A (approx half the atoms)
    lz_full = oriented_unit.cell[2, 2]
    z_cut = lz_full / 2.0
    
    # To avoid cutting *through* a plane, we move z_cut slightly if needed.
    # Or just select atoms < z_cut.
    grain_a_atoms = oriented_unit[positions[:, 2] < z_cut]
    
    # 5. Generate Grain B (Mirror of A)
    # Mirror operation: x->x, y->y, z-> -z
    grain_b_atoms = grain_a_atoms.copy()
    pos_b = grain_b_atoms.get_positions()
    pos_b[:, 2] = -pos_b[:, 2]
    grain_b_atoms.set_positions(pos_b)
    
    # 6. Stack A and B
    # We need to place B on top of A with correct spacing.
    # Find Top of A and Bottom of B (currently B is flipped, so its "bottom" is -z_max_A)
    
    z_a = grain_a_atoms.positions[:, 2]
    z_b = grain_b_atoms.positions[:, 2]
    
    z_max_a = np.max(z_a)
    z_min_b = np.min(z_b) # This is -z_max_a roughly
    
    # We want z_min_new_B = z_max_a + d_spacing
    # Ideally d_spacing is the bulk spacing.
    # We use the median dz from the bulk analysis.
    d_target = d_spacing_measured
    
    shift_z = z_max_a + d_target - z_min_b
    grain_b_atoms.translate([0, 0, shift_z])
    
    # 7. Determine Final Cell Height
    # We need the periodic boundary (Top B -> Bottom A) to also satisfy spacing.
    # z_max_b_new = np.max(grain_b_atoms.positions[:, 2])
    # z_min_a = np.min(z_a)
    # Target Lz = z_max_b_new - z_min_a + d_target
    
    z_max_b_new = np.max(grain_b_atoms.positions[:, 2])
    z_min_a = np.min(z_a)
    
    Lz_new = z_max_b_new - z_min_a + d_target
    
    # 8. Create Final Supercell
    final_cell = oriented_unit.cell.copy()
    final_cell[2, 2] = Lz_new
    
    # Combine
    twin_supercell = grain_a_atoms + grain_b_atoms
    twin_supercell.set_cell(final_cell)
    twin_supercell.set_pbc(True)
    
    # 9. Cleanup & Overlap Removal
    # Since we constructed it with explicit spacing, overlaps should only occur
    # if atoms are laterally coincident (identical X,Y) and Z is flipped identically.
    # Or if d_target was too small.
    # We check for clashes.
    
    twin_supercell.wrap()
    
    dists = get_distances(twin_supercell.get_positions(), cell=twin_supercell.cell, pbc=True)
    dist_mat = dists[1]
    np.fill_diagonal(dist_mat, 100.0)
    
    # Ti-Ti bond ~ 2.9 A.
    # 1.5 A is a safe hard-sphere clash limit.
    clash_limit = 1.5 
    
    to_remove = set()
    rows, cols = np.where(dist_mat < clash_limit)
    
    for r, c in zip(rows, cols):
        if r < c:
            to_remove.add(c)
            
    final_atoms = Atoms(cell=twin_supercell.cell, pbc=True)
    for i, atom in enumerate(twin_supercell):
        if i not in to_remove:
            final_atoms.append(atom)
            
    log_lines.append(f"Construction Details:")
    log_lines.append(f"  Grain A atoms: {len(grain_a_atoms)}")
    log_lines.append(f"  Grain B atoms: {len(grain_b_atoms)}")
    log_lines.append(f"  Z-shift applied: {shift_z:.4f} A")
    log_lines.append(f"  Final Lz: {Lz_new:.4f} A")
    log_lines.append(f"  Overlaps removed: {len(to_remove)}")
    log_lines.append(f"  Final Atom Count: {len(final_atoms)}")
    
    # 10. Write Output
    out_base = "Ti_5811_compact_3d"
    write(os.path.join(output_dir, f"{out_base}.cif"), final_atoms)
    write(os.path.join(output_dir, f"{out_base}.extxyz"), final_atoms)
    
    with open(log_file, 'w') as f:
        f.write("\n".join(log_lines))
        
    print(f"Generated {out_base}. Log: {log_file}")

if __name__ == "__main__":
    build_compact_twin()