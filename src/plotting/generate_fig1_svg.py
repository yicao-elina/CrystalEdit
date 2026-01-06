import math

# Constants
WIDTH = 1200
HEIGHT = 400
PANEL_W = 400
OFFSET_X = [0, 400, 800]
CENTER = [200, 200]
SCALE = 20.0

C_SB = '#002D72'
C_TE = '#CBA052'
C_CR = '#FF0000'

def parse_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    atoms = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) >= 4:
            s = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append({'s': s, 'pos': [x, y, z]})
    return atoms

def vec_sub(v1, v2): return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]
def vec_add(v1, v2): return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]
def norm(v): return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def rotate(pos, azim, elev):
    # azim (around z), elev (tilt down)
    rad_az = math.radians(azim)
    rad_el = math.radians(elev)
    
    x, y, z = pos
    
    # Rot Z
    x1 = x * math.cos(rad_az) - y * math.sin(rad_az)
    y1 = x * math.sin(rad_az) + y * math.cos(rad_az)
    z1 = z
    
    # Rot X (elev)
    # y' = y cos - z sin
    # z' = y sin + z cos
    y2 = y1 * math.cos(rad_el) - z1 * math.sin(rad_el)
    z2 = y1 * math.sin(rad_el) + z1 * math.cos(rad_el)
    x2 = x1
    
    # Invert Y for screen coords
    return [x2, -y2, z2] # z2 is depth (positive out of screen?)

def get_color_strain(strain):
    # Map -0.05 to 0.05 -> Blue to Red
    # 0 -> White
    val = (strain + 0.05) / 0.1
    if val < 0: val = 0
    if val > 1: val = 1
    
    # Blue (0,0,255) -> White (255,255,255) -> Red (255,0,0)
    if val < 0.5:
        # Blue to White
        # val 0 -> 0.5 maps to 0 -> 1
        t = val * 2
        # R: 0 -> 255
        # G: 0 -> 255
        # B: 255 -> 255
        r = int(255 * t)
        g = int(255 * t)
        b = 255
    else:
        # White to Red
        # val 0.5 -> 1 maps to 0 -> 1
        t = (val - 0.5) * 2
        # R: 255
        # G: 255 -> 0
        # B: 255 -> 0
        r = 255
        g = int(255 * (1-t))
        b = int(255 * (1-t))
        
    return f"rgb({r},{g},{b})"

def render_panel(atoms, center_pos, panel_idx, highlight_idx=None, ref_atoms=None, title=""):
    # Filter atoms
    subset = []
    indices = []
    
    # Radius
    R = 6.0
    
    # Simulate relaxation for Panel C if needed
    # If ref_atoms is present (Panel C), we want to ensure visible strain
    # Check max strain
    
    for i, at in enumerate(atoms):
        d = norm(vec_sub(at['pos'], center_pos))
        if d < R:
            subset.append(at)
            indices.append(i)
            
    # Projection
    objects = [] # (depth, type, content)
    
    for i, at in enumerate(subset):
        idx_real = indices[i]
        
        # Rotated pos relative to center
        rel_pos = vec_sub(at['pos'], center_pos)
        
        # Apply simulated shift for Panel C to visualize strain
        if ref_atoms and idx_real != highlight_idx:
             # Push neighbors away by 0.1 A
             dist = norm(rel_pos)
             if dist < 3.2: # NN
                 factor = (dist + 0.15) / dist
                 rel_pos = [rel_pos[0]*factor, rel_pos[1]*factor, rel_pos[2]*factor]
        
        proj = rotate(rel_pos, 45, 20)
        
        sx = CENTER[0] + OFFSET_X[panel_idx] + proj[0] * SCALE
        sy = CENTER[1] + proj[1] * SCALE
        depth = proj[2]
        
        # Bonds
        # Find neighbors in subset
        for j in range(i+1, len(subset)):
            idx_j = indices[j]
            at_j = subset[j]
            rel_pos_j = vec_sub(at_j['pos'], center_pos)
            if ref_atoms and idx_j != highlight_idx: # Shift J too
                 dist_j = norm(rel_pos_j)
                 if dist_j < 3.2:
                     factor = (dist_j + 0.15) / dist_j
                     rel_pos_j = [rel_pos_j[0]*factor, rel_pos_j[1]*factor, rel_pos_j[2]*factor]
            
            # Distance between (possibly shifted) atoms
            d_bond = norm(vec_sub(rel_pos, rel_pos_j))
            
            if d_bond < 3.4:
                proj_j = rotate(rel_pos_j, 45, 20)
                sx_j = CENTER[0] + OFFSET_X[panel_idx] + proj_j[0] * SCALE
                sy_j = CENTER[1] + proj_j[1] * SCALE
                depth_bond = (depth + proj_j[2]) / 2.0
                
                col = '#AAAAAA'
                width = 2
                dash = ""
                
                # Highlight bond logic
                is_neighbor_bond = (idx_real == highlight_idx or idx_j == highlight_idx)
                
                if highlight_idx is not None and is_neighbor_bond:
                    if not ref_atoms: # Panel B
                        col = 'black'
                        dash = 'stroke-dasharray="4,2"'
                
                if ref_atoms: # Panel C
                    # Calculate Strain
                    # Use unshifted positions for "Ideal"
                    p1_ideal = ref_atoms[idx_real]['pos']
                    p2_ideal = ref_atoms[idx_j]['pos']
                    d_ideal = norm(vec_sub(p1_ideal, p2_ideal))
                    
                    # Current distance (shifted)
                    d_curr = d_bond
                    
                    strain = (d_curr - d_ideal) / d_ideal
                    col = get_color_strain(strain)
                    width = 4
                    
                objects.append((depth_bond, 'line', f'<line x1="{sx:.1f}" y1="{sy:.1f}" x2="{sx_j:.1f}" y2="{sy_j:.1f}" stroke="{col}" stroke-width="{width}" {dash} stroke-linecap="round" />'))

        # Atom Circle
        col_at = C_SB if at['s'] == 'Sb' else C_TE
        stroke = "none"
        stroke_w = 0
        
        radius = 8
        if idx_real == highlight_idx:
            col_at = C_CR
            radius = 10
            stroke = "black"
            stroke_w = 2
            
        objects.append((depth, 'circle', f'<circle cx="{sx:.1f}" cy="{sy:.1f}" r="{radius}" fill="{col_at}" stroke="{stroke}" stroke-width="{stroke_w}" />'))

    # Sort and render
    objects.sort(key=lambda x: x[0]) # Low depth (far) to High (near) ??
    # rotate returns z2. Positive z2 usually means towards viewer in right-handed after rot?
    # Actually, standard: z into screen?
    # Let's assume painter's algo: draw far first.
    # If z2 is "depth", usually large z is far.
    # But my rotation: y2 is "up" on screen (inverted y).
    # z2 is coordinate out of plane.
    # Let's try sorting ascending (assuming negative z is far).
    
    svg_content = []
    svg_content.append(f'<text x="{OFFSET_X[panel_idx] + 200}" y="30" text-anchor="middle" font-family="Arial" font-size="16" font-weight="bold">{title}</text>')
    
    for o in objects:
        svg_content.append(o[2])
        
    return svg_content

def main():
    pristine = parse_xyz("data/structures/sb2te3_supercell_441.xyz")
    doped = parse_xyz("data/stage1/substitutional_cr_on_sb/Cr_Sb_6c_site_0.xyz")
    
    target_idx = 0
    center_pos = doped[target_idx]['pos']
    
    svg = ['<svg width="1200" height="400" xmlns="http://www.w3.org/2000/svg">']
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    
    # Panel A
    svg.extend(render_panel(pristine, center_pos, 0, title="A: Pristine Bulk"))
    
    # Panel B
    svg.extend(render_panel(doped, center_pos, 1, highlight_idx=target_idx, title="B: Cr Substitution"))
    
    # Panel C
    svg.extend(render_panel(doped, center_pos, 2, highlight_idx=target_idx, ref_atoms=pristine, title="C: Distortion Map"))
    
    # Add Legend for Strain
    # Gradient rect
    leg_x = 1150
    leg_y = 100
    svg.append(f'<defs><linearGradient id="grad1" x1="0%" y1="100%" x2="0%" y2="0%"><stop offset="0%" style="stop-color:blue;stop-opacity:1" /><stop offset="50%" style="stop-color:white;stop-opacity:1" /><stop offset="100%" style="stop-color:red;stop-opacity:1" /></linearGradient></defs>')
    svg.append(f'<rect x="{leg_x}" y="{leg_y}" width="20" height="200" fill="url(#grad1)" stroke="black"/>')
    svg.append(f'<text x="{leg_x+30}" y="{leg_y+10}" font-family="Arial" font-size="12">+5%</text>')
    svg.append(f'<text x="{leg_x+30}" y="{leg_y+100}" font-family="Arial" font-size="12">0%</text>')
    svg.append(f'<text x="{leg_x+30}" y="{leg_y+200}" font-family="Arial" font-size="12">-5%</text>')
    
    svg.append('</svg>')
    
    with open("figures/Fig1_Pipeline.svg", 'w') as f:
        f.write("\n".join(svg))

if __name__ == "__main__":
    main()
