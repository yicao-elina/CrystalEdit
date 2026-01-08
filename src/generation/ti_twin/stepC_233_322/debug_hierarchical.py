import os
import numpy as np
from ase.io import read, write
from ase import Atoms
from ase.build import make_supercell
from ase.geometry import get_distances

# Configuration
input_cif = "data/structures/ti_twin/Ti-bcc.cif"
# We'll output to a debug directory
outdir = "data/structures/ti_twin/stepC_233_332_twins/debug"
os.makedirs(outdir, exist_ok=True)

def remove_local_overlaps(atoms, frac_center, width=0.03, cutoff=1.35):
    scaled = atoms.get_scaled_positions()
    mask = np.abs(scaled[:, 2] - frac_center) < width
    idx = np.where(mask)[0]

    if len(idx) == 0:
        return atoms, 0

    sub = atoms[idx]
    # Check if get_distances works with sub-selection or needs full cell context
    # It calculates distances between positions in 'sub', using 'atoms.get_cell()'
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

def debug_build():
    print("Debug Hierarchical Builder")
    
    # 1. Load bulk
    bulk = read(input_cif)
    
    # 2. Check {233} matrix P
    # Rows of P are new basis vectors in terms of old basis vectors
    P233 = np.array([
        [ 2,  3, -3],   # A1?
        [ 3, -2, -3],   # A2?
        [ 2,  3,  3]    # A3 (Normal)?
    ])
    # Check orthogonality
    v1 = P233[0]
    v2 = P233[1]
    v3 = P233[2]
    
    print(f"v1.v2 = {np.dot(v1,v2)}") # 6 -6 +9 != 0? Wait.
    # [2, 3, -3] . [3, -2, -3] = 6 - 6 + 9 = 9. Not orthogonal.
    # [3, -2, -3] . [2, 3, 3] = 6 - 6 - 9 = -9. Not orthogonal.
    # [2, 3, -3] . [2, 3, 3] = 4 + 9 - 9 = 4. Not orthogonal.
    
    # The matrix P provided in the script is NOT an orthogonal basis for a cubic box.
    # It defines a triclinic cell.
    # Twin construction by simple Z-mirroring relies on Z being normal to the twin plane
    # AND X,Y being in the twin plane.
    # If the basis is not orthogonal, z-mirroring (1.0 - z) does NOT correspond to a reflection across the plane 
    # unless the c-axis is perpendicular to a, b.
    
    oriented = make_supercell(bulk, P233)
    cell = oriented.get_cell()
    print("Cell angles:", cell.cellpar()[3:])
    
    # If angles are not 90, simple fractional mirroring is invalid for creating a coherent twin boundary.
    # Correct P matrix needed.
    
    # {332} check (Secondary)
    # The script just mirrors a subsection of the {233} twin.
    # It claims to make a {332} twin.
    # Mirroring a {233} crystal across a plane parallel to {233} creates a {233} twin.
    # To create a {332} secondary twin, one must define the secondary twin plane orientation *relative* to the primary.
    # The secondary twin plane is {332} (or {112} of the twin?).
    # It must be a different plane.
    # The current script mirrors across Z (which is the {233} plane).
    # So it just reverses the twin back to original orientation?
    # That makes a Twin-Matrix-Twin structure, but all boundaries are parallel to {233}.
    # It is NOT a hierarchical twin system with different planes.
    
    print("Conclusion: Current script creates parallel twins, not hierarchical structure with different planes.")
    print("And the basis is likely non-orthogonal.")

if __name__ == "__main__":
    debug_build()
