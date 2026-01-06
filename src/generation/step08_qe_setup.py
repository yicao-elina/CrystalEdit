import os

def parse_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Lattice
    lattice = []
    comment = lines[1]
    if 'Lattice="' in comment:
        lat_str = comment.split('Lattice="')[1].split('"')[0]
        v = [float(x) for x in lat_str.split()]
        lattice = [v[0:3], v[3:6], v[6:9]]
    
    atoms = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) >= 4:
            atoms.append({'s': parts[0], 'pos': [float(x) for x in parts[1:4]]})
    return atoms, lattice

def generate_qe_input(filename, atoms, lattice, task_name):
    # Standard QE Template for Sb2Te3
    # PBE functional, ONCV pseudos assumed
    species = sorted(list(set(a['s'] for a in atoms)))
    
    with open(filename, 'w') as f:
        f.write(f"&CONTROL\n")
        f.write(f"    calculation = 'relax',\n")
        f.write(f"    prefix = '{task_name}',\n")
        f.write(f"    outdir = './out/',\n")
        f.write(f"    pseudo_dir = '../pseudo/',\n")
        f.write(f"    forc_conv_thr = 1.0d-4,\n")
        f.write(f"    nstep = 100,\n")
        f.write(f"/\n")
        
        f.write(f"&SYSTEM\n")
        f.write(f"    ibrav = 0,\n")
        f.write(f"    nat = {len(atoms)},")
        f.write(f"    ntyp = {len(species)},")
        f.write(f"    ecutwfc = 50.0,\n")
        f.write(f"    ecutrho = 200.0,\n")
        f.write(f"    occupations = 'smearing',\n")
        f.write(f"    smearing = 'mv',\n")
        f.write(f"    degauss = 0.01,\n")
        # Add magnetism if Cr is present
        if 'Cr' in species:
            f.write(f"    nspin = 2,\n")
            # Starting magmoms
            for i, s in enumerate(species):
                mag = 3.0 if s == 'Cr' else 0.0
                f.write(f"    starting_magnetization({i+1}) = {mag},\n")
        f.write(f"/\n")
        
        f.write(f"&ELECTRONS\n")
        f.write(f"    conv_thr = 1.0d-8,\n")
        f.write(f"    mixing_beta = 0.3,\n")
        f.write(f"/\n")
        
        f.write(f"&IONS\n")
        f.write(f"    ion_dynamics = 'bfgs',\n")
        f.write(f"/\n")
        
        f.write(f"CELL_PARAMETERS (angstrom)\n")
        for v in lattice:
            f.write(f"    {v[0]:12.8f} {v[1]:12.8f} {v[2]:12.8f}\n")
            
        f.write(f"ATOMIC_SPECIES\n")
        # Placeholder pseudos
        for s in species:
            mass = 51.996 if s == 'Cr' else (121.76 if s == 'Sb' else 127.60)
            f.write(f"    {s:<2} {mass:7.3f} {s}.upf\n")
            
        f.write(f"ATOMIC_POSITIONS (angstrom)\n")
        for a in atoms:
            f.write(f"    {a['s']:<2} {a['pos'][0]:12.8f} {a['pos'][1]:12.8f} {a['pos'][2]:12.8f}\n")
            
        f.write(f"K_POINTS (automatic)\n")
        # k-grid scaling
        # Supercell is 17x17x31. 
        # k-grid 2x2x1 or 4x4x1? 
        # For 4x4x1 supercell, 2x2x1 k-points is usually enough.
        f.write(f"    2 2 1 0 0 0\n")

def main():
    targets = [
        ("2.1_slab_Te", "data/structures/sb2te3_slab_5QL_Te_term.xyz"),
        ("2.2_fault", "data/structures/sb2te3_faulted_stacking.xyz"),
        ("3.1_twin", "data/structures/sb2te3_twin_boundary.xyz"),
        ("3.2_gb", "data/structures/sb2te3_gb.xyz"),
        ("cr_doped", "data/stage1/substitutional_cr_on_sb/Cr_Sb_6c_site_0.xyz")
    ]
    
    os.makedirs("data/qe_inputs", exist_ok=True)
    
    for name, path in targets:
        if os.path.exists(path):
            atoms, lattice = parse_xyz(path)
            out_file = f"data/qe_inputs/{name}.in"
            generate_qe_input(out_file, atoms, lattice, name)
            print(f"Generated QE input: {out_file}")

if __name__ == "__main__":
    main()
