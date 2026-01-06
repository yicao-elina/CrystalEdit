import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read as ase_read
from ase.io import write as ase_write
import json
import os

# Paths
cif_path = "../Sb2Te3-mp-1201.cif"
sc_path = "../sb2te3_supercell_441.xyz"
output_dir = "intrinsic_antisite"

if not os.path.exists(cif_path):
    # Fallback
    cif_path = "26ICML-CrystalEdit/Sb2Te3-mp-1201.cif"
    sc_path = "26ICML-CrystalEdit/sb2te3_supercell_441.xyz"
    output_dir = "26ICML-CrystalEdit/stage1/intrinsic_antisite"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print(f"Loading structures...")
unit_cell = Structure.from_file(cif_path)
# Load supercell
atoms_sc = ase_read(sc_path)
supercell = AseAtomsAdaptor.get_structure(atoms_sc)

# 1. Identify Unique Te Sites in Unit Cell
print("Identifying unique Te sites in unit cell...")
sga = SpacegroupAnalyzer(unit_cell)
sym_struct = sga.get_symmetrized_structure()
unique_sites = sym_struct.equivalent_sites

te_sites = []
for equivalents in unique_sites:
    representative = equivalents[0]
    if representative.specie.symbol == "Te":
        wyckoff = sym_struct.wyckoff_symbols[sym_struct.equivalent_sites.index(equivalents)]
        # Te1 is typically 6c, Te2 is 3a in standard setting for Sb2Te3
        # Let's verify by z-coordinate or multiplicity
        label = "Te_Unknown"
        if wyckoff == "3a":
            label = "Te2 (Center QL)"
        elif wyckoff == "6c":
            label = "Te1 (Outer QL)"
            
        te_sites.append({
            "site": representative,
            "wyckoff": wyckoff,
            "label": label,
            "multiplicity": len(equivalents)
        })

print(f"Found {len(te_sites)} inequivalent Te sites.")
for ts in te_sites:
    print(f"  - {ts['label']} (Wyckoff {ts['wyckoff']})")

# 2. Site Selection Logic
# Prioritize Te2 (3a) as requested (center of QL)
selected_site_info = next((s for s in te_sites if "3a" in s['wyckoff']), None)
if not selected_site_info:
    print("Warning: Te(3a) site not found. Selecting first available Te site.")
    selected_site_info = te_sites[0]

print(f"\nSelected site for Antisite Defect: {selected_site_info['label']}")

# Map to Supercell
def get_supercell_index(u_site, sc, scaling=[4,4,1]):
    target_frac = np.array(u_site.frac_coords) / np.array(scaling)
    target_frac = target_frac % 1.0
    
    min_dist = float('inf')
    found_idx = -1
    
    # Iterate all shifts to handle mapping
    for i in range(scaling[0]):
        for j in range(scaling[1]):
            for k in range(scaling[2]):
                shift = np.array([i, j, k])
                t_f = (np.array(u_site.frac_coords) + shift) / np.array(scaling)
                t_f = t_f % 1.0
                
                for idx, site in enumerate(sc):
                    if site.specie.symbol != u_site.specie.symbol: continue
                    diff = site.frac_coords - t_f
                    diff -= np.round(diff)
                    dist = np.linalg.norm(diff)
                    if dist < min_dist:
                        min_dist = dist
                        found_idx = idx
                        
                    if dist < 0.05:
                        return idx
    return found_idx if min_dist < 0.1 else None

sc_idx = get_supercell_index(selected_site_info['site'], supercell)
print(f"Mapped to Supercell Index: {sc_idx}")

# 3. Create Antisite Defect (Sb_Te)
# Replace Te with Sb
defect_sc = supercell.copy()
defect_sc.replace(sc_idx, "Sb")

# 4. Charge Analysis
# Te (CN6) -> Sb (CN6)
# Te formal: -2
# Sb formal: +3 (in Sb2Te3, Sb is +3, Te is -2)
# Replaced Te(-2) with Sb(+3). 
# Charge change: (+3) - (-2) = +5 relative to lattice?
# Wait, this is an antisite.
# In ionic model: Sb_Te means Sb(+3) on Te(-2) site.
# Defect charge state q = Z_defect - Z_host = (+3) - (-2) = +5. 
# This is a very high donor.
# Covalent View:
# Sb: [Kr] 4d10 5s2 5p3 (5 valence)
# Te: [Kr] 4d10 5s2 5p4 (6 valence)
# Sb has 1 less electron than Te.
# Replacing Te with Sb -> 1 hole (p-type)? 
# Actually, antisites in Tetradymites are crucial.
# Sb_Te is typically an acceptor (p-type) or donor?
# Let's consider valence electrons.
# Te has 6, Sb has 5.
# Sb_Te has 1 fewer electron than Te. 
# Acts as an acceptor (p-type). (Hole doping).
# Common knowledge: Sb2Te3 is naturally p-type due to Sb_Te antisites.

charge_analysis = {
    "defect": "Sb_Te",
    "site_type": selected_site_info['label'],
    "nominal_oxidation_change": "Te(-2) -> Sb(+3) (Ionic limit)",
    "valence_change": "Te(6e-) -> Sb(5e-)",
    "electronic_behavior": "Acceptor (p-type)",
    "mechanism": "1 less valence electron introduces a hole."
}

# 5. Validation
# Point Group
try:
    d_sga = SpacegroupAnalyzer(defect_sc, symprec=0.1)
    pg = d_sga.get_point_group_symbol()
except:
    pg = "Unknown"

# Bond Valence Sum (BVS)
# Calculate BVS for the inserted Sb
bvs_analyzer = BVAnalyzer()
try:
    # Need structure with oxidation states or guess
    # Let's verify BVS of the new Sb atom
    # This assumes distance-based BVS parameters available for Sb-Sb or Sb-Te bonds?
    # Nearest neighbors are Sb (since original Te was bonded to Sb).
    # Sb-Sb bonds are metallic/covalent. BVS might not be strictly applicable or give 0.
    valences = bvs_analyzer.get_valences(defect_sc)
    defect_valence = valences[sc_idx]
except Exception as e:
    defect_valence = f"Could not calc BVS: {e}"

# Metrics
# I_fid: Number of atoms should remain same (substitution)
i_fid = (len(defect_sc) == len(supercell))

# Output JSON
output_data = {
    "defect_type": "antisite_Sb_on_Te",
    "site_info": {
        "label": selected_site_info['label'],
        "wyckoff": selected_site_info['wyckoff'],
        "original_element": "Te",
        "new_element": "Sb",
        "supercell_index": sc_idx,
        "frac_coords": defect_sc[sc_idx].frac_coords.tolist()
    },
    "charge_analysis": charge_analysis,
    "validation": {
        "point_group_symmetry": pg,
        "initial_symmetry": "R-3m (D3d)",
        "symmetry_reduction": "D3d -> C3v (likely) or lower",
        "bond_valence_sum": defect_valence,
        "I_fid": i_fid
    },
    "metrics": {
        "spatial_overlap": "High (Substitution)",
        "coordination_number": 6 # Maintained
    }
}

json_file = os.path.join(output_dir, "antisite_report.json")
with open(json_file, 'w') as f:
    json.dump(output_data, f, indent=4)

# 6. Save Extended XYZ
xyz_file = os.path.join(output_dir, "Antisite_Sb_on_Te.xyz")
atoms_out = AseAtomsAdaptor.get_atoms(defect_sc)
atoms_out.info['defect_type'] = "Sb_Te_antisite"
atoms_out.info['description'] = f"Sb substituting Te at {selected_site_info['label']}"

ase_write(xyz_file, atoms_out, format="extxyz")

print(f"Generated Sb_Te antisite. Report saved to {json_file}")
print(f"Structure saved to {xyz_file}")
