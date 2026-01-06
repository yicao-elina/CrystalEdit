import os
import numpy as np
from ase.io import read, write
from ase.build import surface
from ase import Atoms
from ase.geometry import get_distances

# Configuration
input_cif = 'data/structures/ti_twin/Ti-bcc.cif'
output_dir = 'data/structures/ti_twin/stepB_5811_twins'
log_file = 'docs/logs/stepB_5811_twins.log'

# Constants
# (5 8 11) is the twin plane.
# In BCC, the twin plane is usually {112}.
# (5 8 11) is a high index plane.
# The user specified (5 8 11). We proceed with that.
miller_indices = (5, 8, 11)

def build_twin(scale_name, target_z, target_xy_min, vacuum=10.0):
    """
    Builds a twin slab.
    scale_name: 'DFT' or 'MD'
    target_z: Thickness of the crystal part (excluding vacuum)
    target_xy_min: Minimum lateral dimension
    """
    
    # 1. Load primitive
    bulk = read(input_cif)
    
    # 2. Create Surface oriented to (5 8 11)
    # ase.build.surface returns a slab with z perpendicular to surface.
    # It finds the minimal periodic unit cell in x,y.
    # Vacuum is added later.
    slab_prim = surface(bulk, miller_indices, layers=1, vacuum=0)
    
    # Analyze dimensions
    cell_par = slab_prim.cell.cellpar()
    a_len, b_len = cell_par[0], cell_par[1]
    thickness_layer = cell_par[2] # This might be small for 1 layer
    
    # ASE surface 1 layer thickness is usually just the span of atoms.
    # We need to stack layers to reach target_z.
    # Determine thickness per layer (d_spacing).
    # d_hkl = a / sqrt(h^2 + k^2 + l^2) for cubic.
    # Check calculated d vs geometry.
    # (5 8 11) -> 25+64+121 = 210.
    d_spacing = bulk.cell.cellpar()[0] / np.sqrt(210)
    
    # We want half the thickness for Grain A, half for Grain B?
    # Or total thickness `target_z`.
    # Let's target `target_z / 2` for the base slab, then mirror.
    target_half_z = target_z / 2.0
    
    # Estimate layers needed
    # ASE 'layers' argument in surface() controls the repetition in Z direction (sort of).
    # But surface() logic for BCC high index is complex.
    # It repeats the primitive basis?
    # Let's just repeat the slab_prim in Z.
    # slab_prim vector c is often not orthogonal if not carefully set, 
    # but surface() usually tries to make it orthogonal or at least valid for a slab.
    # Let's check orthogonality.
    
    # Better approach:
    # Use surface() with enough layers to get half thickness.
    # Or repeat.
    # Let's find out how many 'layers' give approx d_spacing.
    # It's safer to generate a thicker slab directly if we know the count.
    # Number of planes ~ target_half_z / d_spacing.
    n_planes = int(np.ceil(target_half_z / d_spacing))
    
    slab = surface(bulk, miller_indices, layers=n_planes, vacuum=0)
    slab.center(vacuum=0)
    
    # 3. Lateral Expansion
    # We need lateral dims >= target_xy_min
    nx = int(np.ceil(target_xy_min / a_len))
    ny = int(np.ceil(target_xy_min / b_len))
    
    # For DFT, if target_xy_min is small (e.g. 0), we stick to nx=1, ny=1 
    # UNLESS the minimal cell is tiny.
    # From previous knowledge, minimal (5 8 11) is ~40x22 A.
    # If DFT lateral target is not specified, we usually keep minimal periodic.
    # But prompt says "DFT sizing: Slab ~30 A thickness". 
    # It doesn't enforce XY. Step A had 30x30x30.
    # If we assume 30A box is desired, we might not need to expand XY since 40x22 is close to 30^2 area. 
    # We'll stick to 1x1 for DFT if it fits reasonable bounds, or follow the target_xy_min if given.
    # Let's define target_xy_min for DFT as 0 (keep minimal) and MD as 100.
    
    if nx < 1: nx = 1
    if ny < 1: ny = 1
    
    slab = slab.repeat((nx, ny, 1))
    
    # 4. Generate Twin
    # Grain A: The current slab.
    # Grain B: Mirror of Grain A.
    # Mirror plane: Z = max(z of A) + d_spacing/2 ?
    # Standard twin boundary is a reflection.
    # We reflect atoms coordinates.
    
    grain_a = slab.copy()
    pos_a = grain_a.get_positions()
    z_max_a = np.max(pos_a[:, 2])
    
    # Reflection plane height
    # Ideally halfway between the last plane of A and first of B.
    # Separation should be d_spacing?
    # mirror_z = z_max_a + d_spacing / 2.0
    # Actually, let's just reflect Z -> -Z and then shift.
    
    grain_b = grain_a.copy()
    pos_b = grain_b.get_positions()
    pos_b[:, 2] *= -1
    grain_b.set_positions(pos_b)
    
    # Shift B to sit on top of A
    z_min_b = np.min(grain_b.positions[:, 2])
    # Target: z_min_new = z_max_a + separation
    # Separation for twin? usually related to d_spacing.
    # Let's start with predicted d_spacing.
    shift_z = z_max_a + d_spacing - z_min_b
    grain_b.translate([0, 0, shift_z])
    
    # Combine
    structure = grain_a + grain_b
    
    # 5. Remove duplicates
    # "Remove duplicates at interface within tolerance 0.05 A"
    # Identify atoms at the interface.
    # Interface is around z_max_a.
    # We check all pairs.
    
    # Use ASE get_distances to find close atoms
    # Since we built by stacking, duplicates should only be possible if the mirror operation 
    # placed an atom exactly on top of another (e.g. if an atom was on the mirror plane).
    # Check distances.
    # We only care about periodic in X, Y. Z is not periodic yet.
    # Set provisional cell
    structure.set_pbc([True, True, False])
    # Ensure cell covers the Z height
    z_top = np.max(structure.positions[:, 2])
    z_bot = np.min(structure.positions[:, 2])
    h_c = z_top - z_bot + vacuum
    cell = structure.get_cell()
    cell[2, 2] = h_c
    structure.set_cell(cell)
    # Re-center Z
    structure.center(vacuum=vacuum/2.0, axis=2)
    
    # Identify duplicates
    dists = get_distances(structure.get_positions(), cell=structure.get_cell(), pbc=[True, True, False])
    # dists[1] is distance matrix
    dist_mat = dists[1] + np.eye(len(structure)) * 10.0 # avoid self
    
    to_remove = []
    tolerance = 0.05
    rows, cols = np.where(dist_mat < tolerance)
    
    # Set of indices to remove (remove higher index)
    for r, c in zip(rows, cols):
        if r < c:
            to_remove.append(c)
    
    to_remove = sorted(list(set(to_remove)), reverse=True)
    n_removed = len(to_remove)
    if to_remove:
        del structure[to_remove]
    
    # 6. Finalize
    # Ensure vacuum
    structure.center(vacuum=vacuum/2.0, axis=2)
    
    return structure, n_removed, nx, ny, n_planes

def main():
    log_lines = []
    log_lines.append("Step B: (5 8 11) Twin Supercell Generation")
    log_lines.append("="*40)
    
    # DFT
    log_lines.append("\nGenerating DFT Scale:")
    # Target ~30A total thickness (crystals).
    # Minimal XY usually sufficient for DFT surface calcs unless explicitly asked for bulk supercell size.
    # Since inputs were the bulk supercells (30^3), we'll aim for substantial size.
    # But 40x22 (minimal) is already > 30x30 area. So 1x1 is fine.
    try:
        dft_twin, n_rem_dft, nx_dft, ny_dft, n_layers_dft = build_twin('DFT', target_z=30.0, target_xy_min=0.0, vacuum=10.0)
        
        write(os.path.join(output_dir, 'Ti_5811_DFT.cif'), dft_twin)
        write(os.path.join(output_dir, 'Ti_5811_DFT.extxyz'), dft_twin)
        
        log_lines.append(f"  - Repeats XY: {nx_dft}x{ny_dft}")
        log_lines.append(f"  - Layers per grain: {n_layers_dft}")
        log_lines.append(f"  - Duplicates removed: {n_rem_dft}")
        log_lines.append(f"  - Final Atoms: {len(dft_twin)}")
        log_lines.append(f"  - Cell: {dft_twin.cell.cellpar()}")
    except Exception as e:
        log_lines.append(f"  - Error: {e}")
        print(e)

    # MD
    log_lines.append("\nGenerating MD Scale:")
    # Target >100A XY, >100A Z (total crystal).
    try:
        md_twin, n_rem_md, nx_md, ny_md, n_layers_md = build_twin('MD', target_z=100.0, target_xy_min=100.0, vacuum=10.0)
        
        write(os.path.join(output_dir, 'Ti_5811_MD.cif'), md_twin)
        write(os.path.join(output_dir, 'Ti_5811_MD.extxyz'), md_twin)
        
        log_lines.append(f"  - Repeats XY: {nx_md}x{ny_md}")
        log_lines.append(f"  - Layers per grain: {n_layers_md}")
        log_lines.append(f"  - Duplicates removed: {n_rem_md}")
        log_lines.append(f"  - Final Atoms: {len(md_twin)}")
        log_lines.append(f"  - Cell: {md_twin.cell.cellpar()}")
    except Exception as e:
        log_lines.append(f"  - Error: {e}")
        print(e)

    with open(log_file, 'w') as f:
        f.write('\n'.join(log_lines))
    print(f"Done. Logs in {log_file}")

if __name__ == "__main__":
    main()
