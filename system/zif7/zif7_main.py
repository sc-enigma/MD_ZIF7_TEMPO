import sys
import numpy as np
import pickle

sys.path.append('../components')
from atom import Atom, mol2_to_atoms, count_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, compose_itp_files
from pbc import apply_pbc, calculate_lattice_vectors

from zif7_utils import remove_oxygen, define_zif7_atom_types, define_zif7_atom_names, transform_lattice
from zif7_params import get_zif7_params

# STEP 1. Read data
# Set lower bound in Mercury calculate packing = 0.0
a     = 22.989
b     = 22.989
c     = 15.763
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 120.00 / 180 * np.pi

a_np     = 24.60
b_np     = 21.10
c_np     = 15.70
alpha_np = 84.80 / 180 * np.pi
beta_np  = 86.10 / 180 * np.pi
gamma_np = 128.00 / 180 * np.pi

bounds_a, bounds_b, bounds_c = [1.0, 3.0], [1.0, 3.0], [1.0, 3.0]

vec_a, vec_b, vec_c = calculate_lattice_vectors(a, b, c, alpha, beta, gamma)
atoms = mol2_to_atoms(read_mol2_file('__zif7_source.mol2'))
atoms = remove_oxygen(atoms)

# STEP 2. Apply periodic boundary conditions
atoms = apply_pbc(atoms, a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)

# STEP 3. Define atom types and names
atoms = define_zif7_atom_types(atoms)
atoms = define_zif7_atom_names(atoms)

# STEP 4.a. Write .gro and .mol2 files for large pore structure
with open('__tmp/atoms_zif7_lp.pickle', 'wb') as handle:
    pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)
write_gro_file(atoms, 'zif7_lp.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
write_mol2_file(atoms, 'zif7_lp.mol2', a, b, c, alpha, beta, gamma, skip_long_bonds=True)

# STEP 4.b. Write .gro file for nano pore structure
atoms = transform_lattice(atoms,
                          a, b, c, alpha, beta, gamma,
                          a_np, b_np, c_np, alpha_np, beta_np, gamma_np)
with open('__tmp/atoms_zif7_np.pickle', 'wb') as handle:
    pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)
write_gro_file(atoms, 'zif7_np.gro', a_np, b_np, c_np, alpha_np, beta_np, gamma_np, bounds_a, bounds_b, bounds_c)
write_mol2_file(atoms, 'zif7_np.mol2', a_np, b_np, c_np, alpha_np, beta_np, gamma_np, skip_long_bonds=True)

# STEP 5. Write .itp files
mass, charge, bond_params, angle_params, dihedral_params = get_zif7_params()
write_atoms(atoms, charge, mass, 'ZIF', 'atoms.itp')
write_bonds(atoms, bond_params, 'bonds.itp')
write_angles(atoms, angle_params, 'angles.itp')
write_dihedrals(atoms, dihedral_params, 'dihedrals.itp')
compose_itp_files(['atomtypes.itp',
                   'moleculetype.itp', 'atoms.itp', 'bonds.itp', 'angles.itp', 'dihedrals.itp'], 'zif7.itp')







