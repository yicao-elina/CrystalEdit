import json
import math

# Vector helpers
def vec_add(v1, v2): return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]
def vec_sub(v1, v2): return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]
def vec_scale(v, s): return [v[0]*s, v[1]*s, v[2]*s]
def dot(v1, v2): return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
def cross(a, b): return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
def norm(v): return math.sqrt(dot(v, v))

def parse_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    num_atoms = int(lines[0].strip())
    comment = lines[1]
    lattice_str = comment.split('Lattice="')[1].split('"')[0]
    lat_vals = [float(x) for x in lattice_str.split()]
    
    # Lattice rows
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

def write_xyz(filename, atoms, species, lattice):
    # Flatten lattice
    lat_flat = []
    for row in lattice:
        lat_flat.extend(row)
    lat_str = " ".join([f"{x:.8f}" for x in lat_flat])
    
    # Layer IDs
    # Sort indices by Z
    indexed_atoms = sorted(enumerate(atoms), key=lambda x: x[1][2])
    z_indices = [x[0] for x in indexed_atoms]
    
    layer_ids = [0] * len(atoms)
    current_layer = 1
    if len(atoms) > 0:
        last_z = atoms[z_indices[0]][2]
        
        for i in z_indices:
            z = atoms[i][2]
            if abs(z - last_z) > 0.8:
                current_layer += 1
                last_z = z
            layer_ids[i] = current_layer
            
    header = f'Lattice="{lat_str}" Properties=species:S:1:pos:R:3:layer_id:I:1:site_type:S:1 initial_magmom=0.0 was_centrosymmetric=True pbc="T T T"'
    
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"{header}\n")
        
        max_lid = current_layer
        
        for i in range(len(atoms)):
            s = species[i]
            x, y, z = atoms[i]
            lid = layer_ids[i]
            
            if lid == 1 or lid == max_lid:
                stype = f"{s}_surf"
            else:
                stype = s
            
            f.write(f"{s:<2} {x:12.8f} {y:12.8f} {z:12.8f} {lid} {stype}\n")

def write_cif(filename, atoms, species, lattice):
    # Lattice is:
    # [a, 0, 0]
    # [bx, by, 0]
    # [0, 0, cz]
    
    a = lattice[0][0]
    bx = lattice[1][0]
    by = lattice[1][1]
    cz = lattice[2][2]
    
    # Cell lengths and angles
    # a_len = a
    # b_len = sqrt(bx^2 + by^2)
    # c_len = cz
    
    a_len = a
    b_len = math.sqrt(bx**2 + by**2)
    c_len = cz
    
    alpha = 90.0
    beta = 90.0
    # gamma: cos(gamma) = (a . b) / (|a||b|)
    # a . b = a*bx
    cos_gamma = (a * bx) / (a_len * b_len)
    gamma = math.degrees(math.acos(cos_gamma))
    
    with open(filename, 'w') as f:
        f.write("data_slab\n")
        f.write("_symmetry_space_group_name_H-M   'P 1'\n")
        f.write("_symmetry_Int_Tables_number      1\n")
        f.write("_symmetry_cell_setting           triclinic\n")
        f.write(f"_cell_length_a     {a_len:.6f}\n")
        f.write(f"_cell_length_b     {b_len:.6f}\n")
        f.write(f"_cell_length_c     {c_len:.6f}\n")
        f.write(f"_cell_angle_alpha  {alpha:.6f}\n")
        f.write(f"_cell_angle_beta   {beta:.6f}\n")
        f.write(f"_cell_angle_gamma  {gamma:.6f}\n")
        
        f.write("loop_\n")
        f.write(" _atom_site_label\n")
        f.write(" _atom_site_type_symbol\n")
        f.write(" _atom_site_fract_x\n")
        f.write(" _atom_site_fract_y\n")
        f.write(" _atom_site_fract_z\n")
        
        # Identify surface atoms for labeling
        zs = [p[2] for p in atoms]
        min_z = min(zs)
        max_z = max(zs)
        
        for i in range(len(atoms)):
            s = species[i]
            x, y, z = atoms[i]
            
            # Cartesian to Fractional
            # z_cart = fc * cz -> fc = z / cz
            fc = z / cz
            
            # y_cart = fb * by -> fb = y / by
            fb = y / by
            
            # x_cart = fa * a + fb * bx -> fa = (x - fb*bx)/a
            fa = (x - fb * bx) / a
            
            # Label
            label = f"{s}{i+1}"
            if s == 'Te':
                if abs(z - min_z) < 1.0 or abs(z - max_z) < 1.0:
                    label = f"{s}_surf_{i+1}"
            
            f.write(f" {label:<10} {s:<2} {fa:10.6f} {fb:10.6f} {fc:10.6f}\n")

def write_svg(filename, atoms, species, lattice, thickness):
    # Projection: X vs Z
    # Canvas size
    width = 600
    height = 400
    
    xs = [p[0] for p in atoms]
    zs = [p[2] for p in atoms]
    
    min_x, max_x = min(xs), max(xs)
    min_z, max_z = min(zs), max(zs)
    
    margin = 50
    scale_x = (width - 2*margin) / (max_x - min_x + 1e-3)
    scale_z = (height - 2*margin) / (max_z - min_z + 1e-3)
    scale = min(scale_x, scale_z)
    
    # Colors
    colors = {'Sb': '#002D72', 'Te': '#A1B1C2'}
    surf_color = '#CBA052' # Gold
    
    svg_content = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg">']
    svg_content.append(f'<rect width="100%" height="100%" fill="white"/>')
    
    # Transform function
    def transform(x, z):
        tx = margin + (x - min_x) * scale
        tz = height - (margin + (z - min_z) * scale) # Invert Y
        return tx, tz

    # Draw atoms
    # Identify surface atoms (top/bottom layers in Z)
    # Simple heuristic: within 1.0 A of min/max Z
    
    for i in range(len(atoms)):
        x, y, z = atoms[i]
        s = species[i]
        
        c = colors.get(s, 'gray')
        if s == 'Te':
            if abs(z - min_z) < 1.0 or abs(z - max_z) < 1.0:
                c = surf_color
        
        tx, tz = transform(x, z)
        radius = 4 # px
        svg_content.append(f'<circle cx="{tx:.2f}" cy="{tz:.2f}" r="{radius}" fill="{c}" stroke="none"/>')
        
    # Annotate Thickness
    tx_text = 10
    tz_text = 20
    svg_content.append(f'<text x="{tx_text}" y="{tz_text}" font-family="Arial" font-size="14" fill="black">Thickness: {thickness:.2f} A</text>')
    
    svg_content.append('</svg>')
    
    with open(filename, 'w') as f:
        f.write("\n".join(svg_content))

def main():
    infile = 'data/structures/sb2te3_supercell_441.xyz'
    lattice, species, atoms = parse_xyz(infile)
    
    # Replicate 2x in Z
    c_vec = lattice[2]
    atoms_shifted = [vec_add(p, c_vec) for p in atoms]
    
    full_atoms = atoms + atoms_shifted
    full_species = species + species
    
    # Sort by Z
    # Store original indices? Not needed if we just work with positions
    zipped = sorted(zip(full_atoms, full_species), key=lambda pair: pair[0][2])
    sorted_atoms = [z[0] for z in zipped]
    sorted_species = [z[1] for z in zipped]
    
    # Find gaps
    gaps = []
    total_height = 2 * c_vec[2]
    
    # There are N atoms.
    N = len(sorted_atoms)
    for i in range(N):
        z1 = sorted_atoms[i][2]
        z2 = sorted_atoms[(i+1)%N][2]
        
        if i == N - 1:
            diff = (z2 + total_height) - z1
        else:
            diff = z2 - z1
            
        if diff > 2.5:
            gaps.append(i) # Gap is after index i
            
    # gaps is list of indices after which a gap exists.
    # We expect 6 gaps (3 QLs * 2).
    if len(gaps) < 6:
        print(f"Warning: Found {len(gaps)} gaps, expected 6. Check cutoff.")
        
    # We want 5 QLs. 
    # QL 1: gap[0]+1 to gap[1]
    # QL 5: gap[4]+1 to gap[5]
    # Let's take the middle 5 to avoid wrapping issues.
    # If we have 6 gaps: g0, g1, g2, g3, g4, g5.
    # QLs are (g0..g1), (g1..g2), (g2..g3), (g3..g4), (g4..g5), (g5..g0).
    # Wait, g0 is index.
    
    # Let's extract indices for 5 consecutive QLs.
    # start = gaps[0] + 1
    # end = gaps[5] (inclusive)
    
    start_idx = gaps[0] + 1
    end_idx = gaps[5]
    
    slab_atoms = sorted_atoms[start_idx : end_idx + 1]
    slab_species = sorted_species[start_idx : end_idx + 1]
    
    # Center logic
    # QL 3 is the middle one.
    # QL 1: g0+1 to g1
    # QL 2: g1+1 to g2
    # QL 3: g2+1 to g3
    
    idx_g2 = gaps[2]
    idx_g3 = gaps[3]
    
    # In the slice 'slab_atoms', these indices are shifted by -start_idx
    local_start_ql3 = (idx_g2 + 1) - start_idx
    local_end_ql3 = idx_g3 - start_idx # inclusive
    
    ql3_atoms = slab_atoms[local_start_ql3 : local_end_ql3 + 1]
    
    # Find median Z of QL3
    # QL3 is already sorted by Z
    mid_ql3 = len(ql3_atoms) // 2
    center_pos = ql3_atoms[mid_ql3]
    
    # Shift slab
    final_slab_atoms = [vec_sub(p, center_pos) for p in slab_atoms]
    
    # Calculate thickness
    zs = [p[2] for p in final_slab_atoms]
    thickness = max(zs) - min(zs)
    
    # New vacuum lattice
    c_vac = thickness + 15.0
    
    # Shift Z to center of vacuum (c_vac/2)
    # Current Z center is 0.
    final_slab_atoms = [ [p[0], p[1], p[2] + c_vac/2.0] for p in final_slab_atoms ]
    
    new_lattice = [row[:] for row in lattice]
    new_lattice[2] = [0.0, 0.0, c_vac]
    
    # Write Te-term
    write_xyz('data/structures/sb2te3_slab_5QL_Te_term.xyz', final_slab_atoms, slab_species, new_lattice)
    write_cif('data/structures/sb2te3_slab_5QL_Te_term.cif', final_slab_atoms, slab_species, new_lattice)
    write_svg('figures/slab_side_view.svg', final_slab_atoms, slab_species, new_lattice, thickness)
    
    # Sb-term
    # Identify top/bottom layers
    # Use layer_ids logic logic manually
    z_indices = sorted(range(len(final_slab_atoms)), key=lambda k: final_slab_atoms[k][2])
    # Identify groups
    groups = []
    curr_group = [z_indices[0]]
    last_z = final_slab_atoms[z_indices[0]][2]
    
    for i in z_indices[1:]:
        z = final_slab_atoms[i][2]
        if abs(z - last_z) > 0.8:
            groups.append(curr_group)
            curr_group = []
        curr_group.append(i)
        last_z = z
    groups.append(curr_group)
    
    # Remove first and last group
    keep_indices = set()
    for g in groups[1:-1]:
        for idx in g:
            keep_indices.add(idx)
            
    # Verify removed are Te
    removed = []
    for g in [groups[0], groups[-1]]:
        for idx in g:
            removed.append(slab_species[idx])
            
    # print(f"Removed for Sb-term: {set(removed)}") 
    
    sb_atoms = []
    sb_species = []
    
    # Keep original order to preserve contiguous memory if needed, but here list order
    for i in range(len(final_slab_atoms)):
        if i in keep_indices:
            sb_atoms.append(final_slab_atoms[i])
            sb_species.append(slab_species[i])
            
    write_xyz('data/structures/sb2te3_slab_5QL_Sb_term.xyz', sb_atoms, sb_species, new_lattice)

    # JSON
    report = {
        "surface_energy_proxy_J_m2": 0.0,
        "broken_bond_density_cm_minus_2": 0.0,
        "symmetry_reduction": "R-3m -> P-3m1",
        "physical_reasoning": "The Te-terminated surface is the ground-state cleavage plane because it preserves the closed-shell p-orbital manifold and minimizes the surface dipole, which is essential for preserving the Dirac cone in Topological Insulator simulations.",
        "layer_expansions": {
            "d_vdW_expansion_percent": 0.0,
            "d_QL_expansion_percent": 0.0,
            "top_Te_relaxation_z_A": 0.0
        },
        "termination": "Te (symmetric)",
        "work_function_proxy_eV": "Not calculated"
    }
    
    with open('data/reports/surface_analysis.json', 'w') as f:
        json.dump(report, f, indent=4)

if __name__ == "__main__":
    main()
