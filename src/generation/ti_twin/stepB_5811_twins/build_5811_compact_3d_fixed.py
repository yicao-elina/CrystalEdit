import os
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.build import make_supercell
from ase.geometry import get_distances

# =========================
# Configuration
# =========================
input_cif = 'data/structures/ti_twin/Ti-bcc.cif'
output_dir = 'data/structures/ti_twin/stepB_5811_twins'
log_file = 'docs/logs/stepB_5811_twins_compact_fixed.log'

os.makedirs(output_dir, exist_ok=True)
os.makedirs(os.path.dirname(log_file), exist_ok=True)


def build_compact_twin():
    log = []
    log.append("FIXED 3D PERIODIC (5 8 11) TWIN CONSTRUCTION")
    log.append("=" * 60)

    # -------------------------
    # 1. Load bulk
    # -------------------------
    bulk = read(input_cif)
    a0 = bulk.cell.cellpar()[0]
    log.append(f"Lattice constant a0 = {a0:.4f} Å")

    # -------------------------
    # 2. Oriented (5 8 11) cell
    # -------------------------
    P = np.array([
        [1, -2,  1],   # in-plane
        [5,  1, -3],   # in-plane
        [5,  8, 11]    # normal
    ])

    oriented = make_supercell(bulk, P)

    cell_par = oriented.cell.cellpar()
    log.append("Oriented unit cell:")
    log.append(f"  lx = {cell_par[0]:.3f} Å")
    log.append(f"  ly = {cell_par[1]:.3f} Å")
    log.append(f"  lz = {cell_par[2]:.3f} Å")

    # -------------------------
    # 3. Build bulk-sized supercell
    # -------------------------
    supercell = make_supercell(oriented, [[2, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]])

    supercell.set_pbc(True)

    log.append(f"Supercell atoms (single grain): {len(supercell)}")

    # -------------------------
    # 4. FRACTIONAL twin construction
    # -------------------------
    scaled = supercell.get_scaled_positions()

    # Grain A: lower half
    mask_A = scaled[:, 2] < 0.5
    atoms_A = supercell[mask_A]

    # Grain B: twin of A
    atoms_B = atoms_A.copy()
    scaled_B = atoms_B.get_scaled_positions()
    scaled_B[:, 2] = 1.0 - scaled_B[:, 2]   # fractional twin
    atoms_B.set_scaled_positions(scaled_B)

    # Combine
    bicrystal = atoms_A + atoms_B
    bicrystal.set_cell(supercell.get_cell())
    bicrystal.set_pbc(True)

    log.append(f"Atoms after twin duplication: {len(bicrystal)}")

    # -------------------------
    # 5. Local overlap removal (near twin plane only)
    # -------------------------
    scaled_all = bicrystal.get_scaled_positions()
    interface_mask = np.abs(scaled_all[:, 2] - 0.5) < 0.03

    interface_indices = np.where(interface_mask)[0]
    interface_atoms = bicrystal[interface_indices]

    dmat = get_distances(
        interface_atoms.get_positions(),
        cell=bicrystal.get_cell(),
        pbc=True
    )[1]

    np.fill_diagonal(dmat, 10.0)
    cutoff = 1.35  # Å

    remove_local = set()
    rows, cols = np.where(dmat < cutoff)
    for r, c in zip(rows, cols):
        if r < c:
            remove_local.add(interface_indices[c])

    final_atoms = Atoms(cell=bicrystal.get_cell(), pbc=True)
    for i, at in enumerate(bicrystal):
        if i not in remove_local:
            final_atoms.append(at)

    log.append(f"Removed atoms near twin plane: {len(remove_local)}")
    log.append(f"Final atom count: {len(final_atoms)}")

    # -------------------------
    # 6. Sanity checks
    # -------------------------
    scaled_final = final_atoms.get_scaled_positions()
    assert np.all(scaled_final >= -1e-6)
    assert np.all(scaled_final < 1.0 + 1e-6)

    log.append("Sanity check passed: all atoms inside PBC cell")

    # -------------------------
    # 7. Write outputs
    # -------------------------
    out = "Ti_5811_compact_3d_fixed"
    write(os.path.join(output_dir, f"{out}.cif"), final_atoms)
    write(os.path.join(output_dir, f"{out}.extxyz"), final_atoms)

    with open(log_file, "w") as f:
        f.write("\n".join(log))

    print(f"✓ Fixed twin supercell written to {output_dir}")


if __name__ == "__main__":
    build_compact_twin()
