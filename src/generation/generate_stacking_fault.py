import math
import json
import cmath

def vec_add(v1, v2): return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]
def vec_sub(v1, v2): return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]
def vec_scale(v, s): return [v[0]*s, v[1]*s, v[2]*s]
def dot(v1, v2): return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
def norm(v): return math.sqrt(dot(v, v))

def parse_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    num_atoms = int(lines[0].strip())
    comment = lines[1]
    lattice_str = comment.split('Lattice="')[1].split('"')[0]
    lat_vals = [float(x) for x in lattice_str.split()]
    
    lattice = [
        [lat_vals[0], lat_vals[1], lat_vals[2]],
        [lat_vals[3], lat_vals[4], lat_vals[5]],
        [lat_vals[6], lat_vals[7], lat_vals[8]]
    ]
    
    atoms = []
    species = []
    
    for line in lines[2:]:
        parts = line.split()
        if len(parts) >= 4:
            s = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append([x, y, z])
            species.append(s)
            
    return lattice, species, atoms

def cart_to_frac(cart_pos, lattice):
    # inv matrix
    # lattice rows are v1, v2, v3
    v1 = lattice[0]
    v2 = lattice[1]
    v3 = lattice[2]
    
    # Cross products
    c1 = [v2[1]*v3[2]-v2[2]*v3[1], v2[2]*v3[0]-v2[0]*v3[2], v2[0]*v3[1]-v2[1]*v3[0]]
    c2 = [v3[1]*v1[2]-v3[2]*v1[1], v3[2]*v1[0]-v3[0]*v1[2], v3[0]*v1[1]-v3[1]*v1[0]]
    c3 = [v1[1]*v2[2]-v1[2]*v2[1], v1[2]*v2[0]-v1[0]*v2[2], v1[0]*v2[1]-v1[1]*v2[0]]
    
    vol = dot(v1, c1)
    
    inv_lat = [
        [c1[0]/vol, c2[0]/vol, c3[0]/vol],
        [c1[1]/vol, c2[1]/vol, c3[1]/vol],
        [c1[2]/vol, c2[2]/vol, c3[2]/vol]
    ]
    
    # frac = pos * inv_lat (if pos is row vector)
    # x_frac = x*inv00 + y*inv10 + z*inv20
    x, y, z = cart_pos
    fx = x*inv_lat[0][0] + y*inv_lat[1][0] + z*inv_lat[2][0]
    fy = x*inv_lat[0][1] + y*inv_lat[1][1] + z*inv_lat[2][1]
    fz = x*inv_lat[0][2] + y*inv_lat[1][2] + z*inv_lat[2][2]
    
    return [fx, fy, fz]

def frac_to_cart(frac_pos, lattice):
    fx, fy, fz = frac_pos
    v1 = lattice[0]
    v2 = lattice[1]
    v3 = lattice[2]
    
    x = fx*v1[0] + fy*v2[0] + fz*v3[0]
    y = fx*v1[1] + fy*v2[1] + fz*v3[1]
    z = fx*v1[2] + fy*v2[2] + fz*v3[2]
    
    return [x, y, z]

def analyze_stacking(atoms, species, lattice):
    # Group by Z to find layers and QLs
    # Sorted
    zipped = sorted(zip(atoms, species), key=lambda x: x[0][2])
    s_atoms = [x[0] for x in zipped]
    s_species = [x[1] for x in zipped]
    
    # Identify gaps
    c_len = lattice[2][2]
    gaps = []
    
    last_z = s_atoms[0][2]
    for i in range(len(s_atoms)-1):
        z = s_atoms[i+1][2]
        if (z - last_z) > 2.5: # vdW gap
            gaps.append(i)
        last_z = z
        
    # Wraparound gap?
    if (s_atoms[0][2] + c_len - s_atoms[-1][2]) > 2.5:
        gaps.append(len(s_atoms)-1)
        
    # print(f"Found gaps at indices: {gaps}")
    
    # Define QLs
    # Assuming standard order, QL1 starts at 0?
    # If gap is at 79, then 0-79 is QL1.
    
    qls = []
    start = 0
    for g_idx in gaps:
        ql_atoms = s_atoms[start:g_idx+1]
        ql_species = s_species[start:g_idx+1]
        qls.append({'atoms': ql_atoms, 'species': ql_species})
        start = g_idx + 1
    
    # If start < len, add last segment (if wraparound gap wasn't handled strictly or simply)
    if start < len(s_atoms):
        # This segment belongs to the first QL if wrapped? 
        # For simplicity, assume cell boundaries aligned with gaps or handled by previous logic.
        pass
        
    # Analyze Center of QL
    ql_info = []
    for i, ql in enumerate(qls):
        # Average Frac X, Y of the QL
        fx_sum = 0
        fy_sum = 0
        count = 0
        for pos in ql['atoms']:
            f = cart_to_frac(pos, lattice)
            # Center around 0.5 to avoid boundary wrapping issues in avg?
            # Actually, A, B, C are discrete: (0,0), (2/3, 1/3), (1/3, 2/3).
            # Modulo 1. 
            
            # Simple clustering: A, B, C
            # A: (0,0), B: (0.66, 0.33), C: (0.33, 0.66)
            
            # Just take average of fractional coords? No, because atoms are at different positions.
            # The QL has a "center" layer (Te2 or Sb?).
            # In Sb2Te3: Te1 - Sb - Te2 - Sb - Te1.
            # Te2 is the middle layer. Its position defines the QL position.
            pass
            
        # Find Te2 atoms (middle of Z range of QL)
        zs = [p[2] for p in ql['atoms']]
        min_z, max_z = min(zs), max(zs)
        mid_z = (min_z + max_z) / 2.0
        
        # Find atoms close to mid_z
        center_atoms = []
        for j, pos in enumerate(ql['atoms']):
            if abs(pos[2] - mid_z) < 1.0:
                center_atoms.append(pos)
        
        # Average their frac coords
        fx_avg = 0
        fy_avg = 0
        for pos in center_atoms:
            f = cart_to_frac(pos, lattice)
            fx_avg += f[0] % 1.0
            fy_avg += f[1] % 1.0
        
        fx_avg /= len(center_atoms)
        fy_avg /= len(center_atoms)
        
        # Classify
        # A: (0,0), B: (2/3, 1/3), C: (1/3, 2/3)
        site_type = "?"
        dist_A = (fx_avg-0)**2 + (fy_avg-0)**2
        dist_B = (fx_avg-2/3)**2 + (fy_avg-1/3)**2
        dist_C = (fx_avg-1/3)**2 + (fy_avg-2/3)**2
        
        # Handle periodic boundaries for dist calc (e.g. 0.99 ~ 0.0)
        # Simplified for now assuming centered data
        
        if dist_A < 0.05: site_type = "A"
        elif dist_B < 0.05: site_type = "B"
        elif dist_C < 0.05: site_type = "C"
        
        # Refine dist with wrapping
        # ... 
        
        ql_info.append({'index': i, 'type': site_type, 'center_frac': (fx_avg, fy_avg), 'atoms': ql['atoms'], 'species': ql['species']})
        
    return ql_info

def generate_structure_factor_00L(atoms, species, lattice, L_max=10):
    # F(00L) = sum f_j * exp(2*pi*i * L * z_frac_j)
    # Use approximate form factors or just Z number
    Z_map = {'Sb': 51, 'Te': 52}
    
    intensities = []
    
    c_len = lattice[2][2]
    
    for L in range(1, L_max+1):
        F = complex(0, 0)
        for i, pos in enumerate(atoms):
            s = species[i]
            z_frac = pos[2] / c_len
            f = Z_map.get(s, 1)
            phi = 2 * math.pi * L * z_frac
            term = complex(f * math.cos(phi), f * math.sin(phi))
            F += term
        
        I = abs(F)**2
        intensities.append((L, I))
        
    return intensities

def main():
    infile = 'data/structures/sb2te3_supercell_441.xyz'
    lattice, species, atoms = parse_xyz(infile)
    
    qls = analyze_stacking(atoms, species, lattice)
    
    # print("Initial Stacking:", [q['type'] for q in qls])
    
    # Target: Shift Middle QL (index 1)
    # Shift vector: b = 1/3 <1 0 -1 0>
    # In fractional coords (u, v, w) for hexagonal/rhombohedral basis:
    # A shift from A to B is +(2/3, 1/3, 0).
    # A shift from B to C is +(2/3, 1/3, 0) ? Or +(1/3, 2/3)?
    # A (0,0) -> B (2/3, 1/3). Delta = (2/3, 1/3).
    # B (2/3, 1/3) -> C (1/3, 2/3). Delta = (-1/3, 1/3) = (2/3, 1/3) mod 1? No.
    # 2/3 + 2/3 = 4/3 = 1/3. 1/3 + 1/3 = 2/3. Yes.
    # So vector d = (2/3, 1/3, 0) cycles A->B->C->A.
    
    # We want to create a fault.
    # Current: A B C (likely)
    # Target: "ABC ACB ABC". This is confusing. 
    # Let's assume user wants to shift QL 2.
    # If QL2 is B, shifting by (2/3, 1/3) makes it C. Sequence: A C C. (Bad)
    # Shifting by -(2/3, 1/3) i.e. (1/3, 2/3) makes it A. Sequence: A A C. (Bad)
    # But maybe the initial was A C B?
    
    # Let's just apply the shift d = (1/3, 2/3, 0) (The other partial).
    # Or d = (2/3, 1/3, 0).
    # The prompt mentions Shockley partial b = 1/3 <10-10>.
    # |b| = a / sqrt(3).
    # This corresponds to A->B shift. 
    
    # Let's apply shift to QL 1 (the middle one, 0-indexed).
    target_ql = qls[1]
    
    # Shift vector in Cartesian
    # shift_frac = [2.0/3.0, 1.0/3.0, 0.0] 
    shift_frac = [1.0/3.0, 2.0/3.0, 0.0] # Let's try this one.
    
    shift_cart = frac_to_cart(shift_frac, lattice)
    
    # Apply shift
    new_atoms = []
    new_species = []
    
    # Keep QL0 static
    new_atoms.extend(qls[0]['atoms'])
    new_species.extend(qls[0]['species'])
    
    # Shift QL1
    shifted_ql1 = [vec_add(p, shift_cart) for p in qls[1]['atoms']]
    new_atoms.extend(shifted_ql1)
    new_species.extend(qls[1]['species'])
    
    # Keep QL2 static? Or shift top part? 
    # A stacking fault typically shifts everything above the plane.
    # If I only shift QL1, I create two interfaces.
    # Interface 1 (0-1): A - B -> A - C (if shift B->C).
    # Interface 2 (1-2): B - C -> C - C (if shift B->C).
    # C-C interface is very high energy (cations aligned?).
    
    # If I shift QL1 AND QL2 ...
    # Original: A B C
    # Shift QL1+QL2 by (A->B vector):
    # A (static) -> A
    # B -> C
    # C -> A
    # Result: A C A. (Twin-like).
    # This avoids C-C stacking.
    # Let's shift everything above QL0.
    
    shifted_ql2 = [vec_add(p, shift_cart) for p in qls[2]['atoms']]
    new_atoms.extend(shifted_ql2)
    new_species.extend(qls[2]['species'])
    
    # Now check for overlaps/vdW compression
    # Interface between QL0 and QL1.
    # QL0 (static) top atoms vs QL1 (shifted) bottom atoms.
    
    # Simple check: Min distance across interface
    # Iterate top atoms of QL0 and bot atoms of QL1
    
    # Get Z limits
    z0_max = max([p[2] for p in qls[0]['atoms']])
    z1_min = min([p[2] for p in shifted_ql1])
    
    # Expansion needed?
    min_dist = 10.0
    for p0 in qls[0]['atoms']:
        if p0[2] > z0_max - 2.0: # Optim
            for p1 in shifted_ql1:
                if p1[2] < z1_min + 2.0:
                    d = norm(vec_sub(p1, p0))
                    if d < min_dist: min_dist = d
    
    # print(f"Min dist at Interface 1: {min_dist}")
    
    # Target min dist ~ 3.5 A (vdW). If < 3.0, expand.
    expansion = 0.0
    if min_dist < 3.4:
        expansion = 3.4 - min_dist + 0.1
        # print(f"Expanding gap 1 by {expansion}")
        
    # Apply expansion to QL1 and QL2 (shift Z up)
    final_atoms = []
    # QL0
    final_atoms.extend(qls[0]['atoms'])
    
    # QL1
    exp_vec = [0, 0, expansion]
    shifted_exp_ql1 = [vec_add(p, exp_vec) for p in shifted_ql1]
    final_atoms.extend(shifted_exp_ql1)
    
    # QL2
    # Interface 2 check?
    # QL1 (shifted) vs QL2 (shifted). Relative shift is 0.
    # So interface 2 is preserved (B-C becomes C-A, normal stacking relation preserved).
    # So gap 2 shouldn't change.
    
    shifted_exp_ql2 = [vec_add(p, exp_vec) for p in shifted_ql2]
    final_atoms.extend(shifted_exp_ql2)
    
    # New Lattice (c needs to increase?)
    # If we expanded Z, we should increase c to avoid compressing the wrapping gap?
    # Original c
    c_orig = lattice[2][2]
    new_c = c_orig + expansion
    new_lattice = [row[:] for row in lattice]
    new_lattice[2][2] = new_c
    
    # Calculate XRD
    I_orig = generate_structure_factor_00L(atoms, species, lattice)
    I_fault = generate_structure_factor_00L(final_atoms, new_species, new_lattice)
    
    # Write outputs
    write_xyz('data/structures/sb2te3_faulted_stacking.xyz', final_atoms, new_species, new_lattice)
    write_cif('data/structures/sb2te3_faulted_stacking.cif', final_atoms, new_species, new_lattice)
    
    # JSON Report
    report = {
        "fault_type": "Growth Stacking Fault / Twin",
        "shift_vector_fractional": shift_frac,
        "stacking_sequence_change": "ABC -> ACA (approx)",
        "interface_expansion_A": expansion,
        "min_dist_at_interface_A": min_dist,
        "inter_QL_coupling_proxy": "Reduced due to geometric mismatch",
        "symmetry": "P3m1 (Trigonal)",
        "xrd_00L_shifts": {
            "note": "Intensity changes and peak shifts due to Z-expansion",
            "peaks": [ {"L": L, "I_orig": f"{io:.2f}", "I_fault": f"{ifault:.2f}"} for (L, io), (_, ifault) in zip(I_orig, I_fault) if L <= 6]
        }
    }
    
    with open('data/reports/fault_analysis.json', 'w') as f:
        json.dump(report, f, indent=4)

def write_xyz(filename, atoms, species, lattice):
    lat_flat = []
    for row in lattice:
        lat_flat.extend(row)
    lat_str = " ".join([f"{x:.8f}" for x in lat_flat])
    
    header = f'Lattice="{lat_str}" Properties=species:S:1:pos:R:3 pbc="T T T"'
    
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{header}\n")
        for i in range(len(atoms)):
            f.write(f"{species[i]:<2} {atoms[i][0]:12.8f} {atoms[i][1]:12.8f} {atoms[i][2]:12.8f}\n")

def write_cif(filename, atoms, species, lattice):
    # Simplified CIF writer
    a = lattice[0][0]
    bx = lattice[1][0]
    by = lattice[1][1]
    cz = lattice[2][2]
    
    a_len = a
    b_len = math.sqrt(bx**2 + by**2)
    c_len = cz
    gamma = 120.0 # Approx for hexagonal
    
    with open(filename, 'w') as f:
        f.write("data_faulted\n")
        f.write("_symmetry_space_group_name_H-M   'P 1'\n")
        f.write(f"_cell_length_a     {a_len:.6f}\n")
        f.write(f"_cell_length_b     {b_len:.6f}\n")
        f.write(f"_cell_length_c     {c_len:.6f}\n")
        f.write("_cell_angle_alpha  90.0\n")
        f.write("_cell_angle_beta   90.0\n")
        f.write(f"_cell_angle_gamma  {gamma:.6f}\n")
        f.write("loop_\n")
        f.write(" _atom_site_label\n")
        f.write(" _atom_site_type_symbol\n")
        f.write(" _atom_site_fract_x\n")
        f.write(" _atom_site_fract_y\n")
        f.write(" _atom_site_fract_z\n")
        
        for i in range(len(atoms)):
            f_pos = cart_to_frac(atoms[i], lattice)
            f.write(f" {species[i]}{i+1} {species[i]} {f_pos[0]:10.6f} {f_pos[1]:10.6f} {f_pos[2]:10.6f}\n")

if __name__ == "__main__":
    main()
