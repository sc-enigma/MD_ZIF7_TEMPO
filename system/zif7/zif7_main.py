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
from zif7_params_lp import get_zif7_params_lp
from zif7_params_np import get_zif7_params_np

# STEP 1. Read data
# Set lower bound in Mercury calculate packing = 0.0
a     = 22.989
b     = 22.989
c     = 15.763
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 120.00 / 180 * np.pi

a_lp     = 23.948
b_lp     = 21.354
c_lp     = 16.349
alpha_lp = 90.28 / 180 * np.pi
beta_lp  = 93.28 / 180 * np.pi
gamma_lp = 108.41 / 180 * np.pi

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
write_gro_file(atoms, 'zif7_np.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
write_mol2_file(atoms, 'zif7_np.mol2', a, b, c, alpha, beta, gamma, skip_long_bonds=False)

# STEP 4.b. Write .gro file for nano pore structure
atoms = transform_lattice(atoms,
                          a, b, c, alpha, beta, gamma,
                          a_lp, b_lp, c_lp, alpha_lp, beta_lp, gamma_lp)
with open('__tmp/atoms_zif7_np.pickle', 'wb') as handle:
    pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)
write_gro_file(atoms, 'zif7_lp.gro', a_lp, b_lp, c_lp, alpha_lp, beta_lp, gamma_lp, bounds_a, bounds_b, bounds_c)
write_mol2_file(atoms, 'zif7_lp.mol2', a_lp, b_lp, c_lp, alpha_lp, beta_lp, gamma_lp, skip_long_bonds=False)

# STEP 5. Write .itp files
mass, charge, bond_params, angle_params, dihedral_params = get_zif7_params_np()
skip_long_bonds = False
write_atoms(atoms, charge, mass, 'ZIF', 'atoms.itp')
write_bonds(atoms, bond_params, 'bonds_np.itp',             skip_long_bonds)
write_angles(atoms, angle_params, 'angles_np.itp',          skip_long_bonds)
write_dihedrals(atoms, dihedral_params, 'dihedrals_np.itp', skip_long_bonds)
compose_itp_files(['atomtypes.itp',
                   'moleculetype.itp', 'atoms.itp', 'bonds_np.itp', 'angles_np.itp', 'dihedrals_np.itp'], 'zif7_np.itp')

mass, charge, bond_params, angle_params, dihedral_params = get_zif7_params_lp()
skip_long_bonds = False
write_atoms(atoms, charge, mass, 'ZIF', 'atoms.itp')
write_bonds(atoms, bond_params, 'bonds_lp.itp',             skip_long_bonds)
write_angles(atoms, angle_params, 'angles_lp.itp',          skip_long_bonds)
write_dihedrals(atoms, dihedral_params, 'dihedrals_lp.itp', skip_long_bonds)
compose_itp_files(['atomtypes.itp',
                   'moleculetype.itp', 'atoms.itp', 'bonds_lp.itp', 'angles_lp.itp', 'dihedrals_lp.itp'], 'zif7_lp.itp')






