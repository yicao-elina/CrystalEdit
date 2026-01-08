def convert_xyz(infile, outfile):
    with open(infile, 'r') as fin:
        lines = fin.readlines()
    
    n_atoms = int(lines[0].strip())
    atom_lines = lines[2:]
    
    with open(outfile, 'w') as fout:
        fout.write(f"{n_atoms},\n")
        for line in atom_lines:
            parts = line.split()
            s, x, y, z = parts[0], parts[1], parts[2], parts[3]
            fout.write(f'"{s}", {x}, {y}, {z},\n')

if __name__ == "__main__":
    convert_xyz("data/structures/ti_twin/stepC_233_332_twins/Ti_233_332_hierarchical_fixed.xyz", 
                "src/plotting/atoms.xyz")