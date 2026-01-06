# CrystalEdit: Defect Engineering Framework

This repository contains tools for engineering and analyzing crystallographic defects in $Sb_2Te_3$ and related systems.

## Project Structure

The project is organized as follows:

- **`src/`**: Source code for defect generation, analysis, and plotting.
  - `generation/`: Scripts for generating defect structures (slabs, faults, twins, grain boundaries).
    - `experimental/`: Experimental slab generation scripts.
    - `ti_twin/`: Specific scripts for Titanium twin boundary study.
  - `analysis/`: Scripts for analyzing generated structures (voids, symmetry, validation).
  - `plotting/`: Scripts for generating figures and plots.

- **`data/`**: Generated data and intermediate files.
  - `structures/`: XYZ and CIF structure files.
  - `reports/`: JSON analysis reports.
  - `tables/`: CSV and LaTeX benchmark tables.
  - `qe_inputs/`: Quantum ESPRESSO input files.
  - `stage1/`: Data from Stage 1 (point defects).
  - `stage2/`: Data from Stage 2 (surface slabs).

- **`figures/`**: Generated plots and schematic diagrams (PDF, PNG, SVG).

- **`docs/`**: Documentation and logs.

## Usage

All scripts are designed to be run from the **project root directory**.

### Examples

**Generate a pure surface slab:**
```bash
python3 src/generation/generate_slab_pure.py
```

**Analyze voids in the supercell:**
```bash
python3 src/analysis/analyze_voids_supercell.py
```

**Generate benchmark plots:**
```bash
python3 src/plotting/step09_generate_all_defect_plots.py
```

## Dependencies

- Python 3.x
- NumPy
- ASE (Atomic Simulation Environment)
- Pymatgen
- Matplotlib

Install dependencies via:
```bash
pip install -r requirements.txt
```
