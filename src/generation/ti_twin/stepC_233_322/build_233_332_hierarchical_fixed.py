import os
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.build import make_supercell
from ase.geometry import get_distances

# =============================
# Configuration
# =============================
input_cif = "data/structures/ti_twin/Ti-bcc.cif"
outdir = "data/structures/ti_twin/stepC_233_332_twins"
logfile = "docs/logs/stepC_233_332_fixed.log"

os.makedirs(outdir, exist_ok=True)
os.makedirs(os.path.dirname(logfile), exist_ok=True)

def build_hierarchical_twin():
    log = []
    log.append("FIXED Hierarchical {233} -> {332} Twin Builder")
    log.append("=" * 60)

    # 1. Load Bulk
    bulk = read(input_cif)
    
    # 2. Define Orthogonal Basis for {233}
    # Normal: [2 3 3]
    # In-plane 1: [0 1 -1]
    # In-plane 2: [-3 1 1]
    
    P_ortho = np.array([
        [0, 1, -1],
        [-3, 1, 1],
        [2, 3, 3]
    ])
    
    # Check orthogonality
    # r0.r1 = 0+1-1=0.
    # r0.r2 = 0+3-3=0.
    # r1.r2 = -6+3+3=0.
    # Correct.
    
    oriented = make_supercell(bulk, P_ortho)
    log.append("Established Orthogonal Basis for {233}:")
    log.append(f"  X: [0 1 -1] (Length {np.linalg.norm(oriented.cell[0]):.2f} A)")
    log.append(f"  Y: [-3 1 1] (Length {np.linalg.norm(oriented.cell[1]):.2f} A)")
    log.append(f"  Z: [2 3 3]  (Length {np.linalg.norm(oriented.cell[2]):.2f} A)")
    
    # 3. Create Supercell (Primary A-B)
    # We want a reasonable size.
    # X ~ 4.6 A. Repeat 4x -> 18.4 A.
    # Y ~ 10.8 A. Repeat 2x -> 21.6 A.
    # Z ~ 15.3 A. Repeat 3x -> 45.9 A.
    
    supercell = make_supercell(oriented, [[4,0,0], [0,2,0], [0,3,1]]) # 3 in Z for length?
    # Wait, make_supercell diagonal is safer.
    supercell = make_supercell(oriented, [[4,0,0], [0,2,0], [0,0,3]])
    supercell.set_pbc(True)
    
    log.append(f"Supercell Size: {supercell.cell.cellpar()[:3]}")
    log.append(f"Total Atoms (Single Crystal): {len(supercell)}")
    
    # 4. Construct Primary Twin (A-B)
    # Mirror top half (Z > 0.5)
    
    # Wrap first
    supercell.wrap()
    scaled = supercell.get_scaled_positions()
    z_cut = 0.5
    
    grain_a_mask = scaled[:, 2] < z_cut
    grain_a = supercell[grain_a_mask]
    
    # Grain B: Mirror of A
    # Since basis is orthogonal, we can mirror fractional Z: z -> 1 - z
    # Note: 1-z mirrors across z=0.5.
    
    grain_b = grain_a.copy()
    sb = grain_b.get_scaled_positions()
    sb[:, 2] = 1.0 - sb[:, 2]
    grain_b.set_scaled_positions(sb)
    
    # Correct Z-shift?
    # Mirroring A (0..0.5) via 1-z gives (0.5..1.0).
    # This creates boundaries at 0.5 and 1.0(0.0).
    # We need to check interplanar spacing consistency.
    # We assume 'make_supercell' along Z direction creates periodic planes.
    # Mirroring preserves the spacing at the boundaries if the cut is exactly between planes.
    # Let's verify overlap later.
    
    primary_twin = grain_a + grain_b
    primary_twin.set_cell(supercell.get_cell())
    primary_twin.set_pbc(True)
    
    # Clean overlaps for Primary Twin
    primary_twin = remove_overlaps(primary_twin, log, "Primary Twin")
    
    log.append(f"Primary Twin A-B constructed. Atoms: {len(primary_twin)}")
    
    # 5. Construct Secondary Twin (C inside B)
    # Target: Secondary Twin {332}.
    # We select a domain inside Grain B.
    # Grain B is roughly Z in [0.5, 1.0].
    # Center of C: (0.5, 0.5, 0.75).
    
    # Identify B atoms
    sp = primary_twin.get_scaled_positions()
    # Mask for B domain (approx)
    mask_b = (sp[:, 2] > 0.55) & (sp[:, 2] < 0.95)
    indices_b = np.where(mask_b)[0]
    
    # Define Secondary Plane Normal [3 3 2] in Cartesian
    # We need to transform Miller [3 3 2] to Cartesian.
    # For cubic, [hkl] is direction normal to (hkl).
    # Vector v_sec = [3, 3, 2] * a_lattice?
    # Or just [3, 3, 2] in fractional of the BULK?
    # Yes, direction in bulk.
    # We need to express [3, 3, 2] in the Supercell Cartesian frame.
    
    # Bulk vectors (Cartesian aligned)
    # a1 = [a, 0, 0], a2 = [0, a, 0], a3 = [0, 0, a]
    # So v_sec_cart = [3a, 3a, 2a].
    # Normalize it.
    
    # However, Grain B is TWINNED.
    # Its orientation is reflected relative to bulk.
    # So the local lattice of B is different.
    # If we want a {332} twin *of the B lattice*, we must apply the operation corresponding to {332} in B.
    # Since B is a mirror of A (Bulk), the {332} vector in B is the mirror of the {332} vector in A?
    # Mirror across {233}.
    # Let v_332 = [3, 3, 2].
    # Let n_233 = [2, 3, 3].
    # v_332_in_B = v_332 - 2 * (v_332 . n_233_hat) * n_233_hat.
    # But wait, we just want to create a twin C from B.
    # So we apply the reflection along the {332} plane of B.
    # The normal of that plane in Cartesian space is the vector we just calculated (v_332_in_B).
    
    a0 = bulk.cell.cellpar()[0]
    n_primary = np.array([2, 3, 3])
    n_primary = n_primary / np.linalg.norm(n_primary)
    
    v_sec_A = np.array([3, 3, 2])
    v_sec_A = v_sec_A / np.linalg.norm(v_sec_A)
    
    # Reflect v_sec_A across n_primary to get v_sec_B (the normal in B's frame relative to global)
    # v_reflected = v - 2(v.n)n
    v_sec_B = v_sec_A - 2 * np.dot(v_sec_A, n_primary) * n_primary
    
    # Now we have the normal vector v_sec_B in Cartesian coordinates.
    # We define a Domain C inside B (e.g. a sphere or slab).
    # Let's define a Lenticular domain defined by two planes or just a radius.
    # Sphere for simplicity/stability in debug.
    
    center_C = np.dot(primary_twin.cell, [0.5, 0.5, 0.75]) # Cartesian center
    radius_C = 8.0 # Angstroms
    
    pos = primary_twin.get_positions()
    dists_from_C = np.linalg.norm(pos - center_C, axis=1)
    
    mask_C = (dists_from_C < radius_C) & mask_b
    indices_C = np.where(mask_C)[0]
    
    # Apply Twin Operation to C
    # Reflection across plane passing through center_C with normal v_sec_B.
    # r' = r - 2((r-c).n)n
    
    atoms_C = primary_twin[indices_C]
    pos_C = atoms_C.get_positions()
    
    # Vector from center
    vecs = pos_C - center_C
    
    # Projections
    projs = np.dot(vecs, v_sec_B) # Shape (N,)
    
    # New positions
    # pos_new = pos - 2 * proj * n
    # We need to broadcast n
    offsets = 2 * projs[:, np.newaxis] * v_sec_B
    pos_C_new = pos_C - offsets
    
    # Update positions in the main atoms object
    # Wait, we should update them in place?
    # Yes.
    primary_twin.positions[indices_C] = pos_C_new
    
    # Clean overlaps again (locally around C)
    # The boundary of C will have clashes.
    final_hierarchical = remove_overlaps(primary_twin, log, "Secondary Twin C")
    
    log.append(f"Secondary Twin C inserted (Spherical, r={radius_C} A).")
    log.append(f"Secondary Normal used: {v_sec_B} (Cartesian)")
    log.append(f"Final Atom Count: {len(final_hierarchical)}")
    
    # 6. Save
    write(os.path.join(outdir, "Ti_233_332_hierarchical_fixed.cif"), final_hierarchical)
    write(os.path.join(outdir, "Ti_233_332_hierarchical_fixed.xyz"), final_hierarchical)
    
    with open(logfile, 'w') as f:
        f.write("\n".join(log))
        
    print(f"Generated Ti_233_332_hierarchical_fixed. See {logfile}")

def remove_overlaps(atoms, log, stage_name):
    # Wrap
    atoms.wrap()
    
    dists = get_distances(atoms.get_positions(), cell=atoms.get_cell(), pbc=True)[1]
    np.fill_diagonal(dists, 10.0)
    
    cutoff = 1.4 # Ti
    
    remove = set()
    rows, cols = np.where(dists < cutoff)
    
    for r, c in zip(rows, cols):
        if r < c:
            remove.add(c)
            
    if remove:
        new_atoms = Atoms(cell=atoms.get_cell(), pbc=True)
        for i, at in enumerate(atoms):
            if i not in remove:
                new_atoms.append(at)
        log.append(f"[{stage_name}] Removed {len(remove)} overlapping atoms.")
        return new_atoms
    else:
        log.append(f"[{stage_name}] No overlaps found.")
        return atoms

if __name__ == "__main__":
    build_hierarchical_twin()
