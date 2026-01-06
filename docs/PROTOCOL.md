# PROTOCOL: Stage 2 & 3 Defect Engineering in $Sb_2Te_3$

## 1. Objective
To systematically engineer high-fidelity crystallographic defects in the $Sb_2Te_3$ ($R\bar{3}m$) system, including surface slabs, stacking faults, coherent twin boundaries, and symmetric tilt grain boundaries, for use in downstream AI-integrated materials simulations.

## 2. Mathematical Methods & Logic

### 2.1 Surface Slab Generation
- **Method:** Periodic boundary breakage along the [001] direction.
- **Logic:** Identification of the van der Waals (vdW) gap via neighbor-distance analysis. Cleavage at the $Te^{(1)}-Te^{(1)}$ interface to preserve the $D_{3d}$ point group of individual Quintuple Layers (QLs). 
- **Centering:** Slab centered at $z=0.5$ in fractional coordinates to maintain inversion symmetry and eliminate macroscopic dipoles.

### 2.2 Stacking Fault Induction
- **Transformation:** Applied a Shockley partial dislocation vector $\mathbf{b} = \frac{1}{3} \langle 10\bar{1}0 \rangle$ to the stacking sequence.
- **Logic:** Modified the relative lateral shift of QLs to transition from the pristine **ABCABC** sequence to a faulted **ABCACB** (2H-like) sequence. Rigid-body Z-relaxation was performed to prevent steric clashes.

### 2.3 $\Sigma 3 (10\bar{1}1)$ Twin Boundary
- **Construction:** Domain reflection across the $(10\bar{1}1)$ twinning plane.
- **Logic:** Rotation of the secondary domain to satisfy the Coincidence Site Lattice (CSL) $\Sigma 3$ condition. 
- **Refinement:** Geometric overlap removal with a $1.5$ Å threshold to maintain local stoichiometry and coordination motifs.

### 2.4 $\Sigma 7$ Symmetric Tilt Grain Boundary
- **Geometry:** Rotation around the $[0001]$ tilt axis.
- **CSL Theory:** Identification of the $\Sigma 7$ periodic approximant ($38.21^\circ$) as the optimal match for the target $36.87^\circ$ ($\Sigma 5$ derivative) with $< 1\%$ misfit strain. 
- **Assembly:** Symmetric grain rotation ($ \pm \theta/2$) and overlap removal using a covalent-radius-based density check ($d_{min} < 0.7 \times \sum r_{cov}$). 

## 3. Software & Environment
- **Environment:** Pure Python 3.x (Standard Library).
- **Backend Implementation:** Custom vector and matrix algebra modules for coordinate transformations.
- **Visualization:** SVG-to-PDF conversion using ImageMagick (`convert`).

## 4. Scripts Log
- `step01_slab_gen.py`: (Internal name `generate_slab_pure.py`) Surface engineering.
- `step02_fault_gen.py`: (Internal name `generate_stacking_fault.py`) Planar defect induction.
- `step03_twin_gen.py`: (Internal name `generate_twin_pure.py`) Twin boundary construction.
- `step04_gb_gen.py`: (Internal name `generate_gb_pure.py`) Grain boundary engineering.
- `step05_validation.py`: (Internal name `generate_validation_report.py`) Structural integrity checks.
- `step06_table_gen.py`: (Internal name `generate_table1.py`) Benchmark table creation.
- `step07_fig1_gen.py`: (Internal name `step07_fig1_gen_matplotlib.py`) Multi-panel pipeline schematic using Matplotlib and ASE.
- `step07_fig1_gen_fallback.py`: (Internal name `generate_fig1_svg.py`) SVG-based fallback.
- `step08_qe_setup.py`: Quantum ESPRESSO input generation.
- `step09_all_plots.py`: (Internal name `step09_generate_all_defect_plots.py`) Automated generation of pipeline schematic PDFs for all defect types.
- `step10_stage1_plots.py`: (Internal name `step10_generate_stage1_plots.py`) Batch generation of pipeline schematic PDFs for all Stage 1 defects (interstitials, antisites, vacancies).
- `step11_fig2_gen.py`: (Internal name `step11_generate_fig2.py`) Generation of Figure 2: Statistical benchmark of defect stability vs. spatial diversity (Pareto analysis).

## 5. Stage 4: DFT Simulation Setup
To enable rigorous energetic validation, the generated structures were processed into high-fidelity Quantum ESPRESSO (`pw.x`) input files.
- **Calculation Type:** `relax` (Geometry Optimization).
- **Exchange-Correlation:** GGA-PBE.
- **Pseudopotentials:** Norm-conserving (ONCV) format.
- **Energy Cutoffs:** $E_{cut}^{wfc} = 100$ Ry, $E_{cut}^{rho} = 400$ Ry.
- **K-Point Sampling:** $2 \times 2 \times 1$ Monkhorst-Pack grid (scaled for supercell dimensions to ensure convergence $< 1$ meV/atom).
- **Magnetism:** Spin-polarized ($nspin=2$) with $3.0\ \mu_B$ starting magnetization for Cr-doped systems.
- **Output:** Stored in `26ICML-CrystalEdit/qe_inputs/`.
