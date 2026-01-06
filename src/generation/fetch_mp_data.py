import json
# try:
from mp_api.client import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
# except ImportError as e:
#     print(f"Error importing libraries: {e}")
#     exit(1)

api_key = "SOhpM3MmErxSV7hgObzRRF3HZPjOJaOZ"
material_id = "mp-1201"

print(f"Connecting to Materials Project with ID: {material_id}...")

try:
    with MPRester(api_key) as mpr:
        structure = mpr.get_structure_by_material_id(material_id)
    
    print("Structure retrieved.")
    
    # Convert to ASE
    atoms = AseAtomsAdaptor.get_atoms(structure)
    
    # Write extended xyz
    write("sb2te3.xyz", atoms, format="extxyz")
    print("Written sb2te3.xyz")
    
    # Prepare JSON data
    cell = atoms.get_cell()
    cell_lengths = cell.lengths()
    cell_angles = cell.angles()
    scaled_positions = atoms.get_scaled_positions()
    
    atomic_sites = []
    for i, atom in enumerate(atoms):
        atomic_sites.append({
            "symbol": atom.symbol,
            "position": atom.position.tolist(),
            "scaled_position": scaled_positions[i].tolist()
        })
        
    data = {
        "lattice_parameters": {
            "a": cell_lengths[0],
            "b": cell_lengths[1],
            "c": cell_lengths[2],
            "alpha": cell_angles[0],
            "beta": cell_angles[1],
            "gamma": cell_angles[2]
        },
        "atomic_sites": atomic_sites
    }
    
    with open("structure.json", "w") as f:
        json.dump(data, f, indent=4)
    print("Written structure.json")
    
    print("Done.")

except Exception as e:
    print(f"An error occurred: {e}")
    import traceback
    traceback.print_exc()
