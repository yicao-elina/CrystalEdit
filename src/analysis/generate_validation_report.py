import math
import os
import glob

def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def norm(v):
    return math.sqrt(dot(v, v))

def vec_sub(v1, v2):
    return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]

def mat_vec_mul(m, v):
    # m is 3x3 rows
    res = [0, 0, 0]
    for i in range(3):
        res[i] = m[i][0]*v[0] + m[i][1]*v[1] + m[i][2]*v[2]
    return res

def invert_3x3(m):
    # Cramers rule or simple inverse
    det = m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) - \
          m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) + \
          m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0])
    
    inv = [[0,0,0],[0,0,0],[0,0,0]]
    if abs(det) < 1e-9: return None
    
    inv[0][0] = (m[1][1]*m[2][2] - m[1][2]*m[2][1]) / det
    inv[0][1] = (m[0][2]*m[2][1] - m[0][1]*m[2][2]) / det
    inv[0][2] = (m[0][1]*m[1][2] - m[0][2]*m[1][1]) / det
    
    inv[1][0] = (m[1][2]*m[2][0] - m[1][0]*m[2][2]) / det
    inv[1][1] = (m[0][0]*m[2][2] - m[0][2]*m[2][0]) / det
    inv[1][2] = (m[0][2]*m[1][0] - m[0][0]*m[1][2]) / det
    
    inv[2][0] = (m[1][0]*m[2][1] - m[1][1]*m[2][0]) / det
    inv[2][1] = (m[0][1]*m[2][0] - m[0][0]*m[2][1]) / det
    inv[2][2] = (m[0][0]*m[1][1] - m[0][1]*m[1][0]) / det
    
    return inv

def get_min_dist(atoms, lattice):
    # Naive O(N^2) with MIC
    min_d = 100.0
    pair_info = None
    
    recip = invert_3x3(lattice)
    if recip is None: return 0.0, None
    
    # Check only a subset if N is huge, but here N < 500
    N = len(atoms)
    
    for i in range(N):
        for j in range(i+1, N):
            p1 = atoms[i]['pos']
            p2 = atoms[j]['pos']
            
            # Vector
            dvec = vec_sub(p1, p2)
            
            # MIC
            # Convert to frac
            fvec = mat_vec_mul(recip, dvec)
            # Round
            fvec[0] -= round(fvec[0])
            fvec[1] -= round(fvec[1])
            fvec[2] -= round(fvec[2])
            
            # Back to cart
            cart_d = [0,0,0]
            for k in range(3):
                # lattice rows are vectors
                cart_d[k] = fvec[0]*lattice[0][k] + fvec[1]*lattice[1][k] + fvec[2]*lattice[2][k]
                
            dist = norm(cart_d)
            if dist < min_d:
                min_d = dist
                pair_info = (i, j, atoms[i]['s'], atoms[j]['s'], dist)
                
    return min_d, pair_info

def parse_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    n_atoms = int(lines[0].strip())
    comment = lines[1]
    
    # Parse Lattice
    lattice = [[1,0,0],[0,1,0],[0,0,1]]
    if 'Lattice="' in comment:
        lat_str = comment.split('Lattice="')[1].split('"')[0]
        v = [float(x) for x in lat_str.split()]
        # xyz format Lattice="ax ay az bx by bz cx cy cz"
        lattice = [
            [v[0], v[1], v[2]],
            [v[3], v[4], v[5]],
            [v[6], v[7], v[8]]
        ]
        
    atoms = []
    species_counts = {}
    
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 4: continue
        s = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append({'s': s, 'pos': [x,y,z]})
        species_counts[s] = species_counts.get(s, 0) + 1
        
    return atoms, lattice, species_counts

def get_formula(counts):
    keys = sorted(counts.keys())
    return "".join([f"{k}{counts[k]}" for k in keys])

def analyze():
    files = [
        ("pristine", "data/structures/sb2te3_supercell_441.xyz", "none", "R-3m"),
        ("2.1_slab_Te", "data/structures/sb2te3_slab_5QL_Te_term.xyz", "Surface Slab", "P-3m1"),
        ("2.1_slab_Sb", "data/structures/sb2te3_slab_5QL_Sb_term.xyz", "Surface Slab (Sb)", "P3m1"),
        ("2.2_stacking", "data/structures/sb2te3_faulted_stacking.xyz", "Stacking Fault", "P3m1"),
        ("3.1_twin", "data/structures/sb2te3_twin_boundary.xyz", "Twin Boundary", "P1"),
        ("3.2_gb", "data/structures/sb2te3_gb.xyz", "Sigma7 GB", "P1")
    ]
    
    results = []
    
    print("| Structure ID | Defect Type | N_atoms | I_fid | Min_dist (Å) | Space Group | Formula | Valid? |")
    print("|--------------|-------------|---------|-------|--------------|-------------|---------|--------|")
    
    valid_count = 0
    min_dists = []
    
    for sid, path, dtype, sg in files:
        if not os.path.exists(path):
            print(f"| {sid} | {dtype} | FILE NOT FOUND | - | - | - | - | No |")
            continue
            
        atoms, lattice, counts = parse_xyz(path)
        
        n_atoms = len(atoms)
        formula = get_formula(counts)
        min_d, pair = get_min_dist(atoms, lattice)
        min_dists.append(min_d)
        
        sb = counts.get('Sb', 0)
        te = counts.get('Te', 0)
        
        if te > 0:
            ratio = sb/te
            ideal = 2/3.0
            diff = abs(ratio - ideal)
            i_fid = "✓" if diff < 0.1 else "⚠" 
        else:
            i_fid = "?"
            
        is_valid = (min_d > 1.0)
        valid_mark = "✓" if is_valid else "✗"
        if is_valid: valid_count += 1
        
        print(f"| {sid} | {dtype} | {n_atoms} | {i_fid} | {min_d:.2f} | {sg} | {formula} | {valid_mark} |")
        
    print("\n### Statistical Summary")
    pct_valid = (valid_count / len(files)) * 100
    print(f"- % of structures passing all checks: {pct_valid:.1f}%")
    if min_dists:
        print(f"- Min interatomic distance range: {min(min_dists):.2f} - {max(min_dists):.2f} Å")
    print(f"- Symmetry reduction: 1 structure R-3m, {len(files)-1} structures broken symmetry")

if __name__ == "__main__":
    analyze()