from ase.io import read
import numpy as np
import os

# Path to the CIF file for the supercell
supercell_cif = 'D:/Downloads/CsGeX3/Transition_state/CsGeI3/Random_C_CONTCAR/R_to_state2_Rinv/supercell_to_unitcells/TS.cif'

# Check if the CIF file exists
if not os.path.exists(supercell_cif):
    raise FileNotFoundError(f"The file {supercell_cif} does not exist. Please check the path.")

# Read the supercell from CIF file, specifying the format
supercell = read(supercell_cif, format='cif')

# Ensure the supercell was read correctly
if not hasattr(supercell, 'cell'):
    raise ValueError("The CIF file could not be read correctly as an ASE Atoms object. Check the file format or content.")

# Define the supercell dimensions (e.g., 2x2x2 supercell)
supercell_dimensions = (2, 2, 2)  # adjust as necessary

# Get the unit cell from the supercell's cell parameters
unit_cell_matrix = supercell.cell / supercell_dimensions
print("The unitcell lattice parameters are: ", unit_cell_matrix) 

# Decompose supercell into individual unit cells
tolerance = 10e-5
def decompose_supercell(supercell, dimensions, tolerance):
    unit_cells = []
    inverse_unit_cell_matrix = np.linalg.inv(unit_cell_matrix)
    n_unit_cells = np.prod(dimensions)  # Total number of unit cells
    
    # Create a list to hold fractional coordinates for each unit cell
    for _ in range(n_unit_cells):
        unit_cells.append([])

    # Classify atoms based on their fractional coordinates in the unit cell
    for atom in supercell:
        # Calculate fractional coordinates within the supercell
        fractional_coords = np.dot(atom.position, inverse_unit_cell_matrix)
        
        # Determine the unit cell indices without rounding, using a tolerance for range checks
        cell_index = []
        for coord, dim in zip(fractional_coords, dimensions):
            index = int(coord) if (coord % 1) < (1 - tolerance) else int(coord + 1) % dim
            cell_index.append(index)
        
        # Calculate linear index for the cell
        linear_index = (cell_index[0] * dimensions[1] * dimensions[2] +
                        cell_index[1] * dimensions[2] +
                        cell_index[2])
        
        # Store the fractional coordinates within the unit cell
        unit_cells[linear_index].append(fractional_coords % 1)

    return unit_cells

# Decompose the supercell into unit cells
unit_cells = decompose_supercell(supercell, supercell_dimensions, tolerance)

# Output the fractional coordinates of atoms in each unit cell
for i, unit_cell_atoms in enumerate(unit_cells):
    print(f"\nUnit cell {i + 1} contains atoms at fractional coordinates:")
    for frac_coords in unit_cell_atoms:
        print(frac_coords)
    print(f"Total atoms in unit cell {i + 1}: {len(unit_cell_atoms)}")
