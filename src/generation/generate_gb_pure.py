import math
import json

# ==========================================
# Math Helpers
# ==========================================
def dot(v1, v2):
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

def cross(v1, v2):
    return [
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    ]

def norm(v):
    return math.sqrt(dot(v, v))

def normalize(v):
    l = norm(v)
    if l < 1e-9: return [0,0,0]
    return [v[0]/l, v[1]/l, v[2]/l]

def vec_sub(v1, v2):
    return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]

def vec_add(v1, v2):
    return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]

def vec_scale(v, s):
    return [v[0]*s, v[1]*s, v[2]*s]

def mat_vec_mul(m, v):
    res = [0, 0, 0]
    for i in range(3):
        res[i] = m[i][0]*v[0] + m[i][1]*v[1] + m[i][2]*v[2]
    return res

def rotate_z(pos, theta_rad):
    c = math.cos(theta_rad)
    s = math.sin(theta_rad)
    x, y, z = pos
    nx = x*c - y*s
    ny = x*s + y*c
    return [nx, ny, z]

# ==========================================
# CSL Search
# ==========================================

def get_prim_vectors_hex(a, c):
    v1 = [a, 0.0, 0.0]
    v2 = [-0.5*a, 0.5*a*math.sqrt(3.0), 0.0]
    v3 = [0.0, 0.0, c]
    return [v1, v2, v3]

def find_csl_hex(target_angle_deg, a):
    # Search for integer vectors (n, m) in hex basis
    # V = n*a1 + m*a2
    # Squared length = a^2 * (n^2 + m^2 - nm)  <-- Wait, a1.a2 = a*(-a/2) = -a^2/2? 
    # v1 = (a, 0), v2 = (-a/2, a*rt3/2).
    # v1.v2 = -a^2/2.
    # |V|^2 = n^2 a^2 + m^2 a^2 + 2nm (v1.v2) = a^2 (n^2 + m^2 - nm).
    
    # We want to find a vector V such that the angle between V and its mirror or some other vector matches target.
    # For Symmetric Tilt, we rotate V by theta to V'.
    # In CSL theory for rotation around [001] in Hex:
    # Sigma values are of form n^2 + m^2 + nm? No, that's length squared.
    # Sigma = (n^2 + m^2 - nm) / common_factor?
    # Actually, for rotation theta:
    # tan(theta/2) = ...
    
    # Let's just brute force vectors V1, V2 such that angle(V1, V2) ~ target_angle?
    # No, for STGB, we rotate the crystal.
    # We need a vector V in the lattice that becomes V' after rotation, where V' is ALSO a lattice vector?
    # That defines the CSL.
    # Angle theta is defined by coincidence.
    
    # Common Hex CSL angles:
    # Sigma 7: 21.79 deg, 38.21 deg.
    # Sigma 13: 27.8 deg, 32.2 deg.
    # Sigma 19: 13.17 deg, 46.83 deg.
    
    # Target is 36.87 deg.
    # 38.21 is closest. Difference 1.34 deg.
    # Let's check Sigma 7.
    # Sigma 7 vectors: (3, -1) and (2, 1)?
    # Let's compute angle between (2,1) and (1,2)?
    # |(2,1)|^2 = 4 + 1 - 2 = 3.
    # |(3,1)|^2 = 9 + 1 - 3 = 7.
    
    # Let's search for the "Periodicity Vector" along the boundary.
    # If we rotate by theta, the boundary plane is at theta/2.
    # We need a vector along the boundary in the average lattice.
    # Effectively, we need integer vectors (n,m) in Lattice A and (p,q) in Lattice B such that:
    # V_A (rotated by +theta/2) || V_B (rotated by -theta/2)
    # i.e. V_A || V_B if we undo rotation?
    # This implies V_A and V_B are related by rotation theta.
    # V_B = R(theta) V_A.
    
    # So we search for lattice vectors V_A, V_B such that |V_A| = |V_B| and angle is close to target.
    
    best_fit = None
    min_strain_score = 1000.0 # Combo score
    
    search_range = 5
    
    # Store the absolute best purely by angle difference as a fallback
    fallback_fit = None
    min_angle_diff = 1000.0

    for n1 in range(-search_range, search_range+1):
        for m1 in range(-search_range, search_range+1):
            if n1 == 0 and m1 == 0: continue
            
            # Vector 1
            # Cartesian
            v1x = n1*a - m1*0.5*a
            v1y = m1*0.5*a*math.sqrt(3.0)
            l1 = math.sqrt(v1x**2 + v1y**2)
            
            for n2 in range(-search_range, search_range+1):
                for m2 in range(-search_range, search_range+1):
                    if n2 == 0 and m2 == 0: continue
                    
                    v2x = n2*a - m2*0.5*a
                    v2y = m2*0.5*a*math.sqrt(3.0)
                    l2 = math.sqrt(v2x**2 + v2y**2)
                    
                    # Strain (length mismatch)
                    strain = abs(l1 - l2) / ((l1+l2)/2)
                    # Relaxed tolerance
                    if strain > 0.05: continue
                    
                    # Angle
                    dot_val = v1x*v2x + v1y*v2y
                    cos_theta = dot_val / (l1 * l2)
                    if cos_theta > 1.0: cos_theta = 1.0
                    if cos_theta < -1.0: cos_theta = -1.0
                    angle = math.degrees(math.acos(cos_theta))
                    
                    diff = abs(angle - target_angle_deg)
                    
                    candidate = {
                        'v1': (n1, m1),
                        'v2': (n2, m2),
                        'angle': angle,
                        'l1': l1,
                        'l2': l2,
                        'strain': strain,
                        'cart_v1': [v1x, v1y, 0],
                        'cart_v2': [v2x, v2y, 0]
                    }

                    if diff < min_angle_diff:
                        min_angle_diff = diff
                        fallback_fit = candidate
                    
                    if diff < 10.0: # Within 10 degrees
                        score = diff + strain * 100 # Weight strain heavily
                        if score < min_strain_score:
                            min_strain_score = score
                            best_fit = candidate
                            
    if best_fit is None:
        return fallback_fit
    return best_fit

# ==========================================
# Main
# ==========================================

def main():
    a = 4.30391797
    c = 31.77741530
    
    target_angle = 36.87
    
    csl = find_csl_hex(target_angle, a)
    
    print(f"Nearest CSL found: {csl['angle']:.2f} deg (Target {target_angle})")
    print(f"Indices: {csl['v1']} and {csl['v2']}")
    
    # We will use this angle as the "Stress-free" rotation angle for the bicrystal
    # The actual crystals will be rotated by +/- angle/2.
    
    theta_deg = csl['angle']
    theta_rad = math.radians(theta_deg)
    
    # Periodic vectors along the boundary
    # The vector along the boundary is the mean of V1 and V2?
    # Symmetric tilt: Boundary plane bisects V1 and V2.
    # So if we rotate V1 by +theta/2 and V2 by -theta/2, they align.
    # This aligned vector is the boundary period vector 'L_period'.
    
    # Let's define the Frame:
    # Boundary Plane Normal will be Y axis (Cartesian).
    # Tilt Axis is Z axis (Cartesian) = Lattice c.
    # Boundary Period Vector is X axis (Cartesian).
    
    # Wait, original crystal c is [0,0,1].
    # Tilt axis is [001]. So Z is Tilt Axis.
    # Plane is (hk0). Normal is in XY plane.
    # Period vector is in XY plane.
    
    # Construct Supercell Basis for the Bicrystal
    # L_period length = l1 (approx l2).
    # Since strain < 1%, we average them.
    l_avg = (csl['l1'] + csl['l2']) / 2.0
    
    # We need to build a rectangular box for the Grain.
    # Basis vectors for Grain (before rotation):
    # A_grain = V1 (from CSL search)
    # B_grain = Vector perpendicular to V1 in basal plane?
    # C_grain = Original c axis.
    
    # Find orthogonal vector to V1 in basal plane
    # V1 = (v1x, v1y, 0)
    # V_perp = (-v1y, v1x, 0)
    # Is V_perp a lattice vector?
    # In hex, perpendicular vectors are not always simple integers.
    # We need a periodic supercell.
    # We can multiply V_perp by a factor to make it commensurate.
    # Or just cut a slab and accept non-periodicity in Y?
    # "Periodic Approximant ... periodic along boundary plane".
    # So periodicity in Z (tilt axis) is trivial (c).
    # Periodicity in X (boundary vector) is L_avg.
    # Periodicity in Y (boundary normal) is not required for a slab (we have 2 grains).
    # But for the simulation box to be periodic, we usually stack A-B-A-B.
    # Let's assume a Slab geometry with vacuum in Y.
    
    # Grain A Setup:
    # Rotate such that V1 aligns with Global X.
    # But Grain A is rotated by -theta/2.
    # So V1 (which is at angle alpha in crystal A) needs to end up at angle 0 in global?
    # Wait.
    # Symmetric Tilt: Boundary plane bisects the crystals.
    # If V1 corresponds to (h k l) and V2 to (h k l) symmetry equivalent?
    # In our search, V1=(3,1) and V2=(2,-1)? No.
    # Let's assume we align V1 of A to Global X.
    # And we align V2 of B to Global X.
    # Then we rotate A by 0? And B by 0?
    # But A and B are misoriented by theta.
    # So we rotate A by -theta/2 and B by +theta/2 relative to the "Symmetry Plane".
    
    # Let's keep it simple:
    # 1. Generate supercell of Grain A aligned along V1.
    #    Basis: X' = V1, Y' = Perp, Z' = c.
    # 2. Generate supercell of Grain B aligned along V2.
    #    Basis: X'' = V2, Y'' = Perp, Z'' = c.
    # 3. Rotate Grain A by -theta/2.
    #    Actually, if we defined the supercell such that X is V1, then V1 is along X.
    #    If we rotate by -theta/2, V1 ends up at -theta/2.
    #    We want V1_A and V2_B to meet at the interface?
    #    For STGB, the boundary plane (say XZ) bisects the crystals.
    #    Grain A (left): Lattice rotated such that [hk0] is normal?
    #    Grain B (right): Lattice rotated symmetric.
    #    If the boundary is the interface, the two crystals share a period vector along the boundary.
    #    Let's assume the period vector is the mean of V1 and V2?
    #    No, usually we define the boundary by the lattice plane.
    #    $\Sigma 7$ boundary is often $(1 2 \bar{3} 0)$.
    #    Let's assume the CSL vectors V1 and V2 we found ARE the period vectors along the boundary.
    #    So we want V1_A to align with Y axis (Boundary Direction)?
    #    Let's define:
    #    Boundary Normal = X axis.
    #    Boundary Period = Y axis.
    #    Tilt Axis = Z axis.
    
    #    Grain A: Rotate so V1 aligns with Y. (Rotation 1).
    #             Then Rotate by -theta/2 around Z.
    #    Grain B: Rotate so V2 aligns with Y. (Rotation 2).
    #             Then Rotate by +theta/2 around Z.
    #    Wait, if V1 and V2 are separated by theta in the lattice frame?
    #    Then if we align them both to Y, they are effectively rotated relative to lattice by theta.
    #    If we then rotate A by -theta/2 and B by +theta/2...
    #    Net misorientation is theta.
    #    This seems correct.
    
    # Supercell Construction:
    # We need a box for Grain A that is periodic in Y (length |V1|) and Z (length c).
    # And finite in X (Normal).
    # We need to find a vector perpendicular to V1 in the lattice to define the X-width?
    # Or just cut a shape.
    # For stoichiometry and periodic boundary conditions in Y, we need the box Y-vector to be exactly V1.
    # Box Z-vector is c.
    # Box X-vector? Can be non-periodic if we have vacuum.
    # Let's assume vacuum in X.
    
    # Step 1: Generate Atoms for Grain A in a box defined by V1 and V_perp.
    # V1 = csl['v1'] (indices).
    # Find integer vector V_perp orthogonal to V1?
    # In hex, n1*n2 + m1*m2 - 0.5(...) ?
    # Dot product metric: G = [[a^2, -0.5a^2], [-0.5a^2, a^2]].
    # v1 . v_perp = 0.
    # [n1, m1] G [n2, m2]T = 0.
    # n1(n2 - 0.5m2) + m1(-0.5n2 + m2) = 0
    # n1*n2 - 0.5 n1 m2 - 0.5 m1 n2 + m1 m2 = 0
    # 2 n1 n2 - n1 m2 - m1 n2 + 2 m1 m2 = 0
    # n2 (2 n1 - m1) + m2 (2 m1 - n1) = 0
    # Solution: n2 = -(2 m1 - n1), m2 = (2 n1 - m1).
    # Let's check.
    # (2n1 - m1)(-(2m1-n1)) + (2m1-n1)(2n1-m1) = 0. Yes.
    # So V_perp indices: (-2m1 + n1, 2n1 - m1).
    # Let's calculate this vector for Grain A.
    
    n1, m1 = csl['v1']
    np1 = -2*m1 + n1
    mp1 = 2*n1 - m1
    
    # This V_perp is orthogonal to V1.
    # We can use it to define the width of the grain.
    # Generate atoms in grid defined by V1 and V_perp.
    
    # Grain B similar with V2.
    n2, m2 = csl['v2']
    np2 = -2*m2 + n2
    mp2 = 2*n2 - m2
    
    # Generate points.
    
    # Grain A:
    # Rotate atoms by angle such that V1 aligns with Global Y.
    # Angle of V1 = atan2(v1y, v1x).
    # Rot_A = -phi1 + 90 (to align with Y).
    # Then tilt by -theta/2.
    # Total Rot A = 90 - phi1 - theta/2.
    
    # Grain B:
    # Angle of V2 = atan2(v2y, v2x).
    # Rot_B = -phi2 + 90.
    # Then tilt by +theta/2.
    # Total Rot B = 90 - phi2 + theta/2.
    
    # Verify relative rotation:
    # Rot_B - Rot_A = (-phi2 + phi1) + theta.
    # We know phi2 - phi1 approx theta (or -theta).
    # So (-theta) + theta = 0?
    # Wait.
    # We want final misorientation theta.
    # Current A is aligned Y. Current B is aligned Y.
    # They are "identical" in global frame if we just aligned V1 and V2.
    # But V1 and V2 are different vectors in the lattice.
    # So the lattices are misoriented by the angle between V1 and V2 (which is theta).
    # So if we align V1 to Y and V2 to Y, we have ALREADY rotated the lattices relative to each other by theta.
    # Grain A was rotated by -phi1. Grain B by -phi2.
    # Delta = phi1 - phi2 = -theta.
    # So relative orientation is theta.
    # Do we need the extra +/- theta/2?
    # If we just align V1 to Y and V2 to Y, the boundary plane (XZ? No, we said Normal=X, Period=Y) is the interface.
    # The interface is defined by the line where we joined them.
    # Is this the Symmetric Tilt?
    # Symmetry means the crystal planes (hkl) in A and (hkl) in B are symmetric wrt boundary.
    # If V1 and V2 are symmetry equivalent vectors (e.g. (2,1) and (1,2) in Hex), then yes.
    # In our search, we pick V1, V2 such that |V1|=|V2|.
    # So yes, they are symmetric.
    # So simply aligning V1 to Y and V2 to Y creates the STGB.
    # We don't need additional rotation if V1 and V2 are the symmetric pair.
    
    # Construction:
    # Grain A: Left side (X < 0).
    # Grain B: Right side (X > 0).
    # Join at X=0.
    
    # Coordinates generation:
    # For Grain A:
    #  Basis V1, V_perp, c.
    #  Fill box.
    #  Rotate so V1 -> Y.
    #  Shift so X is negative.
    
    # For Grain B:
    #  Basis V2, V_perp2, c.
    #  Fill box.
    #  Rotate so V2 -> Y.
    #  Shift so X is positive.
    
    unit_atoms = get_atoms_in_cell(a, c)
    
    # Generate A
    atoms_A = generate_grain(a, c, unit_atoms, csl['v1'], (np1, mp1), "A")
    # Rotate A
    phi1 = math.degrees(math.atan2(csl['cart_v1'][1], csl['cart_v1'][0]))
    rot_A = 90 - phi1
    atoms_A = apply_rotation(atoms_A, rot_A)
    
    # Generate B
    atoms_B = generate_grain(a, c, unit_atoms, csl['v2'], (np2, mp2), "B")
    # Rotate B
    phi2 = math.degrees(math.atan2(csl['cart_v2'][1], csl['cart_v2'][0]))
    rot_B = 90 - phi2
    atoms_B = apply_rotation(atoms_B, rot_B)
    
    # Calculate box width (X direction)
    # V_perp length?
    l_perp1 = norm([np1*a - mp1*0.5*a, mp1*0.5*a*math.sqrt(3), 0])
    # We generated e.g. 1 unit of V_perp.
    # The width is l_perp1.
    # Shift A to left. Center at -l_perp1/2.
    # Shift B to right. Center at +l_perp1/2.
    
    # Wait, we need to ensure the CUT is planar.
    # The unit cell generation makes a parallelepiped.
    # When rotated, the "V_perp" side might not be along X?
    # V_perp is orthogonal to V1.
    # Since V1 is along Y, V_perp is along X.
    # So the box is rectangular in XY. Good.
    
    # Shift A
    for at in atoms_A:
        at['pos'][0] -= l_perp1/2.0
    
    # Shift B
    for at in atoms_B:
        at['pos'][0] += l_perp1/2.0
        
    # Combine
    all_atoms = atoms_A + atoms_B
    
    # Overlap Removal
    final_atoms = remove_overlaps(all_atoms)
    
    # Output
    write_output(final_atoms, l_avg, l_perp1*2, c, csl)

def generate_grain(a, c, unit_atoms, v_indices, vp_indices, grain_id):
    # Basis vectors for supercell
    n, m = v_indices
    np, mp = vp_indices
    
    # Cartesian basis
    v_vec = [n*a - m*0.5*a, m*0.5*a*math.sqrt(3), 0]
    vp_vec = [np*a - mp*0.5*a, mp*0.5*a*math.sqrt(3), 0]
    c_vec = [0, 0, c]
    
    # Grid search to fill this basis
    # Since we use integer combinations of prim vectors, we can just iterate prim indices (ni, mi, li)
    # And check if inside the parallelepiped defined by v_vec, vp_vec, c_vec.
    # Range?
    # Max index approx length of v_vec / a.
    
    rng = 4 # safe range
    
    atoms = []
    
    # Inverse basis matrix to check fractional coords in supercell
    # M = [v_vec, vp_vec, c_vec] (cols)
    # We treat Z separately.
    
    det = v_vec[0]*vp_vec[1] - v_vec[1]*vp_vec[0]
    
    prim_a1 = [a, 0, 0]
    prim_a2 = [-0.5*a, 0.5*a*math.sqrt(3), 0]
    
    for ni in range(-rng, rng+1):
        for mi in range(-rng, rng+1):
            # Origin of this unit cell
            org = vec_add(vec_scale(prim_a1, ni), vec_scale(prim_a2, mi))
            
            for sp, frac in unit_atoms:
                # Atom pos
                # pos = org + frac[0]a1 + frac[1]a2 + frac[2]c
                dpos = vec_add(vec_scale(prim_a1, frac[0]), vec_scale(prim_a2, frac[1]))
                pos = vec_add(org, dpos)
                # Z coord
                z = frac[2] * c
                pos[2] = z
                
                # Check bounds in Supercell Frame
                # x, y part
                # x = u * vx + w * vpx
                # y = u * vy + w * vpy
                # Solve for u, w
                # u = (x*vpy - y*vpx) / det
                # w = (y*vx - x*vy) / det
                
                u = (pos[0]*vp_vec[1] - pos[1]*vp_vec[0]) / det
                w = (pos[1]*v_vec[0] - pos[0]*v_vec[1]) / det
                
                if u >= -0.001 and u < 0.999 and w >= -0.001 and w < 0.999:
                    atoms.append({
                        's': sp,
                        'pos': pos,
                        'grain': 1 if grain_id=="A" else 2
                    })
                    
    return atoms

def apply_rotation(atoms, angle_deg):
    rad = math.radians(angle_deg)
    for at in atoms:
        at['pos'] = rotate_z(at['pos'], rad)
    return atoms

def remove_overlaps(atoms):
    # Sort by X to make neighbor search easier? Or just brute force (N is small ~200)
    keep_indices = set(range(len(atoms)))
    
    # Identify boundary: X ~ 0
    # Check atoms with abs(x) < 2.5
    
    candidates = [i for i, a in enumerate(atoms) if abs(a['pos'][0]) < 3.0]
    
    radii = {'Sb': 1.4, 'Te': 1.35} # Covalent radii approx
    
    to_remove = set()
    
    for i in range(len(candidates)):
        idx1 = candidates[i]
        if idx1 in to_remove: continue
        
        for j in range(i+1, len(candidates)):
            idx2 = candidates[j]
            if idx2 in to_remove: continue
            
            a1 = atoms[idx1]
            a2 = atoms[idx2]
            
            # Distance
            d2 = (a1['pos'][0]-a2['pos'][0])**2 + (a1['pos'][1]-a2['pos'][1])**2 + (a1['pos'][2]-a2['pos'][2])**2
            d = math.sqrt(d2)
            
            limit = 0.7 * (radii[a1['s']] + radii[a2['s']])
            
            if d < limit:
                # Remove one. Prefer removing from Grain 2?
                if a1['grain'] == 2:
                    to_remove.add(idx1)
                else:
                    to_remove.add(idx2)
                    
    final = []
    for i in range(len(atoms)):
        if i not in to_remove:
            final.append(atoms[i])
            
    return final

def get_atoms_in_cell(a, c):
    # Same as previous script, standard Wyckoff
    sb_z = 0.3958
    te1_z = 0.2175
    
    atoms = []
    # Te2 (3a)
    atoms.append(('Te', [0,0,0]))
    atoms.append(('Te', [2/3, 1/3, 1/3]))
    atoms.append(('Te', [1/3, 2/3, 2/3]))
    
    # Te1 (6c)
    for shift in [[0,0,0], [2/3,1/3,1/3], [1/3,2/3,2/3]]:
        atoms.append(('Te', vec_add(shift, [0,0,te1_z])))
        atoms.append(('Te', vec_add(shift, [0,0,-te1_z])))
        
    # Sb (6c)
    for shift in [[0,0,0], [2/3,1/3,1/3], [1/3,2/3,2/3]]:
        atoms.append(('Sb', vec_add(shift, [0,0,sb_z])))
        atoms.append(('Sb', vec_add(shift, [0,0,-sb_z])))
        
    # Wrap
    for i in range(len(atoms)):
        p = atoms[i][1]
        atoms[i] = (atoms[i][0], [p[0]%1, p[1]%1, p[2]%1])
        
    return atoms

def write_output(atoms, period_y, width_x, length_z, csl):
    # Lattice: 
    # L1 = [0, period_y, 0] ?? No, usually X is L1.
    # Our Period vector is Y. Width is X.
    # Let's verify output convention.
    # If we want visual Top View in XY, it's fine.
    # L1 = [width_x, 0, 0]
    # L2 = [0, period_y, 0]
    # L3 = [0, 0, length_z]
    
    # Note: width_x includes Vacuum? No, we packed them tight?
    # We generated 1 unit of V_perp for A and B.
    # So width is L_perp + L_perp = 2*L_perp.
    # And we shifted them to touch at 0.
    # We should add vacuum in X to isolate the GB?
    # Or is it a bicrystal with 2 GBs?
    # Periodic in X implies 2 GBs.
    # We constructed it such that outer edges match?
    # Left edge of A (at -L_perp) and Right edge of B (at +L_perp).
    # A is -theta/2. B is +theta/2.
    # At the outer boundary, we wrap from B to A.
    # Difference is theta.
    # So we have a second GB at the boundary.
    # This is standard for periodic bicrystals.
    
    with open("data/structures/sb2te3_gb.xyz", 'w') as f:
        f.write(f"{len(atoms)}\n")
        lat = f"{width_x:.6f} 0.0 0.0 0.0 {period_y:.6f} 0.0 0.0 0.0 {length_z:.6f}"
        header = f'Lattice="{lat}" Properties=species:S:1:pos:R:3:grain_id:I:1:dist_to_GB:R:1 misorientation_angle={csl["angle"]:.2f} axis=[001] Sigma=7 pbc="T T T"'
        f.write(f"{header}\n")
        for a in atoms:
            x, y, z = a['pos']
            dist = abs(x) # Distance to center GB
            # Wrap distance if > width/2?
            f.write(f"{a['s']:<2} {x:12.6f} {y:12.6f} {z:12.6f} {a['grain']} {dist:.4f}\n")
            
    # JSON Report
    report = {
        "Sigma": 7,
        "Target_Angle": 36.87,
        "Actual_Angle": csl['angle'],
        "Misfit_Strain_Percent": csl['strain']*100,
        "Grain_Boundary_Energy_Proxy": "Coordination deficit analysis required",
        "Physical_Reasoning": "The Sigma 7 tilt boundary (nearest to Sigma 5 in hex) creates a periodic array of dislocations. In Sb2Te3, these structural units may alter local Berry curvature.",
        "Stoichiometry": {
            "Sb": len([x for x in atoms if x['s']=='Sb']),
            "Te": len([x for x in atoms if x['s']=='Te'])
        }
    }
    with open("data/reports/gb_report.json", 'w') as f:
        json.dump(report, f, indent=4)
        
    # SVG
    write_svg(atoms, width_x, period_y)

def write_svg(atoms, w, h):
    # Top View (XY)
    # Only plot a slice in Z to avoid clutter?
    # Or project all?
    # Let's plot z in [0, 5] approx (one layer).
    
    slice_atoms = [a for a in atoms if a['pos'][2] < 5.0]
    
    svg = [f'<svg width="800" height="400" xmlns="http://www.w3.org/2000/svg">']
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    
    # Scale
    # min x = -w/2, max x = w/2.
    # min y = -h/2? No y is 0 to h.
    # shift to canvas
    
    sx = 700 / w
    sy = 300 / h
    s = min(sx, sy)
    
    cx = 400
    cy = 200
    
    # Colors
    colors = {'Sb': '#002D72', 'Te': '#CBA052'}
    
    for a in slice_atoms:
        px = cx + a['pos'][0] * s
        py = cy + a['pos'][1] * s # Invert Y?
        
        col = colors.get(a['s'], 'gray')
        stroke = "none"
        width = 0
        if abs(a['pos'][0]) < 4.0:
            stroke = "black"
            width = 2
            
        svg.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="5" fill="{col}" stroke="{stroke}" stroke-width="{width}"/>')
        
    # Annotate
    svg.append(f'<line x1="{cx}" y1="0" x2="{cx}" y2="400" stroke="black" stroke-dasharray="5,5"/>')
    svg.append(f'<text x="{cx+10}" y="20" font-family="Arial">GB Plane</text>')
    
    svg.append('</svg>')
    
    with open("figures/bicrystal_top_view.svg", 'w') as f:
        f.write("\n".join(svg))

if __name__ == "__main__":
    main()
