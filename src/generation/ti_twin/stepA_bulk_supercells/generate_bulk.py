import os
import numpy as np
from ase.io import read, write
from ase.build import make_supercell
import math

# Configuration
input_cif = 'data/structures/ti_twin/Ti-bcc.cif'
output_dir = 'data/structures/ti_twin/stepA_bulk_supercells'
log_file = 'docs/logs/stepA_bulk_supercells.log'

# Targets
target_dft = 30.0 # Angstrom
target_md = 100.0 # Angstrom (10 nm)

def generate_supercell(atoms, target_size):
    """Calculates repeat factors and generates supercell."""
    cell = atoms.get_cell()
    lengths = np.linalg.norm(cell, axis=1)
    
    # Calculate repeats (ceil)
    repeats = np.ceil(target_size / lengths).astype(int)
    
    # Ensure periodic boundaries are maintained by using diagonal matrix if input is aligned
    # or just tuple repeat if using ase.build.make_supercell with simple multiplier
    # atoms.repeat(repeats) works for simple diagonal expansion.
    
    supercell = atoms.repeat(repeats)
    return supercell, repeats

def main():
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Load input
    try:
        prim = read(input_cif)
    except Exception as e:
        print(f"Error reading {input_cif}: {e}")
        return

    log_lines = []
    log_lines.append(f"Step A: Bulk Supercell Generation")
    log_lines.append(f"Input: {input_cif}")
    log_lines.append(f"Primitive Cell: {prim.cell.cellpar()}")
    log_lines.append("-" * 40)

    # 1. DFT Supercell
    dft_super, dft_rep = generate_supercell(prim, target_dft)
    dft_name_base = "bulk_Ti_DFT"
    
    write(os.path.join(output_dir, f"{dft_name_base}.cif"), dft_super)
    write(os.path.join(output_dir, f"{dft_name_base}.extxyz"), dft_super)
    
    log_lines.append(f"DFT Supercell Target: > {target_dft} A")
    log_lines.append(f"Replication Factors: {dft_rep}")
    log_lines.append(f"Final Size: {dft_super.cell.cellpar()[:3]}")
    log_lines.append(f"Number of Atoms: {len(dft_super)}")
    log_lines.append(f"Files: {dft_name_base}.cif, {dft_name_base}.extxyz")
    log_lines.append("-" * 40)

    # 2. MD Supercell
    md_super, md_rep = generate_supercell(prim, target_md)
    md_name_base = "bulk_Ti_MD"
    
    write(os.path.join(output_dir, f"{md_name_base}.cif"), md_super)
    write(os.path.join(output_dir, f"{md_name_base}.extxyz"), md_super)
    
    log_lines.append(f"MD Supercell Target: > {target_md} A")
    log_lines.append(f"Replication Factors: {md_rep}")
    log_lines.append(f"Final Size: {md_super.cell.cellpar()[:3]}")
    log_lines.append(f"Number of Atoms: {len(md_super)}")
    log_lines.append(f"Files: {md_name_base}.cif, {md_name_base}.extxyz")
    log_lines.append("-" * 40)

    # Write log
    with open(log_file, 'w') as f:
        f.write('\n'.join(log_lines))
    
    print(f"Bulk generation complete. Log written to {log_file}")

if __name__ == "__main__":
    main()
