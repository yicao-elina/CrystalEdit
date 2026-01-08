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
logfile = "docs/logs/stepC_233_332.log"

os.makedirs(outdir, exist_ok=True)
os.makedirs(os.path.dirname(logfile), exist_ok=True)

# =============================
# Helper: local overlap removal
# =============================
def remove_local_overlaps(atoms, frac_center, width=0.03, cutoff=1.35):
    scaled = atoms.get_scaled_positions()
    mask = np.abs(scaled[:, 2] - frac_center) < width
    idx = np.where(mask)[0]

    if len(idx) == 0:
        return atoms, 0

    sub = atoms[idx]
    dmat = get_distances(sub.get_positions(),
                          cell=atoms.get_cell(),
                          pbc=True)[1]
    np.fill_diagonal(dmat, 10.0)

    remove = set()
    rows, cols = np.where(dmat < cutoff)
    for r, c in zip(rows, cols):
        if r < c:
            remove.add(idx[c])

    new_atoms = Atoms(cell=atoms.get_cell(), pbc=True)
    for i, at in enumerate(atoms):
        if i not in remove:
            new_atoms.append(at)

    return new_atoms, len(remove)

# =============================
# Main builder
# =============================
def build_233_332_twin():
    log = []
    log.append("Hierarchical {233}->{332} Twin Builder (BCC Ti)")
    log.append("=" * 60)

    # ---- 1. Load bulk ----
    bulk = read(input_cif)
    log.append(f"Bulk atoms: {len(bulk)}")

    # ---- 2. Build {233} oriented cell ----
    P233 = np.array([
        [ 2,  3, -3],   # in-plane
        [ 3, -2, -3],   # in-plane
        [ 2,  3,  3]    # normal {233}
    ])

    cell_233 = make_supercell(bulk, P233)

    # enlarge for bulk-like region
    cell_233 = make_supercell(cell_233,
                              [[2,0,0],
                               [0,2,0],
                               [0,0,1]])
    cell_233.set_pbc(True)

    log.append(f"{'{233}'} oriented atoms: {len(cell_233)}")

    # ---- 3. Primary {233} twin ----
    scaled = cell_233.get_scaled_positions()

    A_mask = scaled[:, 2] < 0.5
    A = cell_233[A_mask]

    B = A.copy()
    scaled_B = B.get_scaled_positions()
    scaled_B[:, 2] = 1.0 - scaled_B[:, 2]
    B.set_scaled_positions(scaled_B)

    twin_233 = A + B
    twin_233.set_cell(cell_233.get_cell())
    twin_233.set_pbc(True)

    twin_233, nrm1 = remove_local_overlaps(twin_233, 0.5)

    log.append(f"Primary twin built, atoms: {len(twin_233)}, removed: {nrm1}")

    # ---- 4. Identify primary twin domain ----
    scaled = twin_233.get_scaled_positions()
    primary_domain_mask = (scaled[:, 2] > 0.25) & (scaled[:, 2] < 0.5)
    primary_domain = twin_233[primary_domain_mask]

    # ---- 5. Secondary {332} twin INSIDE primary twin ----
    # operate in the same fractional frame
    sec = primary_domain.copy()
    scaled_sec = sec.get_scaled_positions()

    # secondary twin plane at center of this domain
    z0 = np.mean(scaled_sec[:, 2])
    scaled_sec[:, 2] = 2*z0 - scaled_sec[:, 2]

    sec.set_scaled_positions(scaled_sec)

    # ---- 6. Merge hierarchical structure ----
    final = Atoms(cell=twin_233.get_cell(), pbc=True)

    primary_indices = set(primary_domain_mask.nonzero()[0])
    for i, at in enumerate(twin_233):
        if i not in primary_indices:
            final.append(at)

    for at in sec:
        final.append(at)

    final, nrm2 = remove_local_overlaps(final, z0)

    log.append(f"Secondary twin inserted, removed: {nrm2}")
    log.append(f"Final atom count: {len(final)}")

    # ---- 7. Sanity check ----
    sf = final.get_scaled_positions()
    assert np.all(sf >= -1e-6) and np.all(sf < 1.0 + 1e-6)

    # ---- 8. Write outputs ----
    write(os.path.join(outdir, "Ti_233_332_hierarchical.cif"), final)
    write(os.path.join(outdir, "Ti_233_332_hierarchical.extxyz"), final)

    with open(logfile, "w") as f:
        f.write("\n".join(log))

    print("✓ Hierarchical {233}->{332} twin generated successfully")

# =============================
if __name__ == "__main__":
    build_233_332_twin()
