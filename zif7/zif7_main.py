import sys
import numpy as np
import pickle

sys.path.append('../components')
from atom import Atom, mol2ToAtoms, countAtoms
from read_utils import readMol2File
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, compose_itp_files
from pbc import apply_PBC, calculate_lattice_vectors

from zif7_utils import remove_oxygen, defineAtomTypes, defineAtomNames, transform_lattice

# Set lower bound in Mercury calculate packing = 0.0
a     = 22.989
b     = 22.989
c     = 15.763
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 120.00 / 180 * np.pi

bounds_a, bounds_b, bounds_c = [1.0, 3.0], [1.0, 3.0], [1.0, 3.0]

vec_a, vec_b, vec_c = calculate_lattice_vectors(a, b, c, alpha, beta, gamma)

# [atom_type] = m
mass = {}
mass['intraFF_Zn'] = 65.360
mass['intraFF_N']  = 14.007
mass['intraFF_C1'] = 12.011
mass['intraFF_C2'] = 12.011
mass['intraFF_C5'] = 12.011
mass['intraFF_C6'] = 12.011
mass['intraFF_H1'] = 1.0080
mass['intraFF_H5'] = 1.0080
mass['intraFF_H6'] = 1.0080

# [atom_type] = q
charge = {}
charge['intraFF_Zn'] = 0.7828
charge['intraFF_N']  =-0.4502
charge['intraFF_C1'] = 0.2432
charge['intraFF_C2'] = 0.1722
charge['intraFF_C5'] =-0.0881
charge['intraFF_C6'] =-0.1982
charge['intraFF_H1'] = 0.0988
charge['intraFF_H5'] = 0.0988
charge['intraFF_H6'] = 0.0988

# [atom_type-atom_type] = [funct, r0, k]
bond_params = {}
bond_params['intraFF_C1-intraFF_H1'] = [1, 0.1090, 154544.408]
bond_params['intraFF_C1-intraFF_N']  = [1, 0.1340, 146958.816]
bond_params['intraFF_C2-intraFF_N']  = [1, 0.1390, 114733.648]
bond_params['intraFF_C2-intraFF_C2'] = [1, 0.1420, 120637.272]
bond_params['intraFF_C2-intraFF_C5'] = [1, 0.1400, 146561.336]
bond_params['intraFF_C5-intraFF_H5'] = [1, 0.1090, 151360.384]
bond_params['intraFF_C5-intraFF_C6'] = [1, 0.1390, 153523.512]
bond_params['intraFF_C6-intraFF_H6'] = [1, 0.1090, 150962.904]
bond_params['intraFF_C6-intraFF_C6'] = [1, 0.1410, 142247.632]
bond_params['intraFF_Zn-intraFF_N']  = [3, 0.1990, 109.0350, 19.8000]

# [atom_type-atom_type-atom_type] = [funct, angle, k]
angle_params = {}
angle_params['intraFF_N-intraFF_C1-intraFF_H1']  = [1, 122.77, 227.610]
angle_params['intraFF_N-intraFF_C1-intraFF_N']   = [1, 114.47, 473.880]
angle_params['intraFF_N-intraFF_C2-intraFF_C5']  = [1, 131.01, 628.939]
angle_params['intraFF_N-intraFF_C2-intraFF_C2']  = [1, 107.68, 510.406]
angle_params['intraFF_C1-intraFF_N-intraFF_C2']  = [1, 105.09, 598.479]
angle_params['intraFF_C2-intraFF_C5-intraFF_H5'] = [1, 121.52, 251.542]
angle_params['intraFF_C2-intraFF_C5-intraFF_C6'] = [1, 117.01, 658.436]
angle_params['intraFF_C5-intraFF_C6-intraFF_C6'] = [1, 121.67, 796.968]
angle_params['intraFF_C5-intraFF_C6-intraFF_H6'] = [1, 119.17, 251.458]
angle_params['intraFF_C6-intraFF_C6-intraFF_H6'] = [1, 119.15, 253.341]
angle_params['intraFF_C6-intraFF_C5-intraFF_H5'] = [1, 121.48, 267.274]
angle_params['intraFF_C2-intraFF_C2-intraFF_C5'] = [1, 120.93, 610.027]
angle_params['intraFF_N-intraFF_Zn-intraFF_N']   = [1, 109.42, 76.1490]
angle_params['intraFF_Zn-intraFF_N-intraFF_C1']  = [1, 126.76, 76.4840]
angle_params['intraFF_Zn-intraFF_N-intraFF_C2']  = [1, 127.94, 78.8270]

# [atom_type-atom_type-atom_type-atom_type] = [funct, angle, k, n]        - periodic
# [atom_type-atom_type-atom_type-atom_type] = [funct, c1, c2, c3, c4, c5] - fourier
dihedral_params = {}
dihedral_params['intraFF_N-intraFF_C2-intraFF_C5-intraFF_C6']  = [9, 180.00, 1.58992, 2]
dihedral_params['intraFF_H1-intraFF_C1-intraFF_N-intraFF_C2']  = [9, 180.00, 9.91608, 2]
dihedral_params['intraFF_C1-intraFF_N-intraFF_C2-intraFF_C2']  = [9, 180.00, 2.88696, 2]
dihedral_params['intraFF_H5-intraFF_C5-intraFF_C2-intraFF_N']  = [9, 180.00, 0.92048, 2]
dihedral_params['intraFF_C2-intraFF_C2-intraFF_C5-intraFF_H5'] = [9, 180.00, 1.08784, 2]
dihedral_params['intraFF_C5-intraFF_C2-intraFF_C2-intraFF_C5'] = [9, 180.00, 2.59408, 2]
dihedral_params['intraFF_H5-intraFF_C5-intraFF_C6-intraFF_C6'] = [9, 180.00, 0.92048, 2]
dihedral_params['intraFF_H6-intraFF_C6-intraFF_C6-intraFF_C5'] = [9, 180.00, 2.34304, 2]
dihedral_params['intraFF_C5-intraFF_C6-intraFF_C6-intraFF_C5'] = [9, 180.00, 4.43504, 2]
dihedral_params['intraFF_N-intraFF_C2-intraFF_C2-intraFF_N']   = [9, 180.00, 4.60240, 2]
dihedral_params['intraFF_N-intraFF_C2-intraFF_C2-intraFF_C5']  = [9, 180.00, 2.05016, 2]
dihedral_params['intraFF_C1-intraFF_N-intraFF_C2-intraFF_C5']  = [9, 180.00, 2.46856, 2]
dihedral_params['intraFF_C2-intraFF_C5-intraFF_C6-intraFF_H6'] = [9, 180.00, 5.18816, 2]
dihedral_params['intraFF_C6-intraFF_C5-intraFF_C2-intraFF_C2'] = [9, 180.00, 0.92048, 2]
dihedral_params['intraFF_C6-intraFF_C6-intraFF_C5-intraFF_C2'] = [9, 180.00, 2.42672, 2]
dihedral_params['intraFF_N-intraFF_C1-intraFF_N-intraFF_C2']   = [9, 180.00, 9.62320, 2]
dihedral_params['intraFF_C1-intraFF_N-intraFF_Zn-intraFF_N']   = [5, 12.63568, 9.45584, 0.33472, 0.00000]
dihedral_params['intraFF_C2-intraFF_N-intraFF_Zn-intraFF_N']   = [5, 14.18376, -8.03328, 0.08368, 0.00000]

# STEP 1. Read data
atoms = mol2ToAtoms(readMol2File('__zif7_source.mol2'))

# STEP 2. Remove oxygen
atoms = remove_oxygen(atoms)

# STEP 3. Apply periodic boundary conditions
atoms = apply_PBC(atoms, a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)

# STEP 4. Define atom types
atoms = defineAtomTypes(atoms)

# STEP 5. Update atom names
atoms = defineAtomNames(atoms)

# STEP 6.a Write .gro and .mol2 files for large pore structure
write_gro_file(atoms, 'zif7_lp.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
write_mol2_file(atoms, 'zif7_lp.mol2', a, b, c, alpha, beta, gamma, skip_long_bonds=True)

# STEP 6.b Write .gro file for nano pore structure
a_np     = 24.60
b_np     = 21.10
c_np     = 15.70
alpha_np = 84.80 / 180 * np.pi
beta_np  = 86.10 / 180 * np.pi
gamma_np = 128.00 / 180 * np.pi

atoms = transform_lattice(atoms,
                          a, b, c, alpha, beta, gamma,
                          a_np, b_np, c_np, alpha_np, beta_np, gamma_np)
write_gro_file(atoms, 'zif7_np.gro', a_np, b_np, c_np, alpha_np, beta_np, gamma_np, bounds_a, bounds_b, bounds_c)
write_mol2_file(atoms, 'zif7_np.mol2', a_np, b_np, c_np, alpha_np, beta_np, gamma_np, skip_long_bonds=True)

# STEP 7. Write .itp files
write_atoms(atoms, charge, mass, 'ZIF', 'atoms.itp')
write_bonds(atoms, bond_params, 'bonds.itp')
write_angles(atoms, angle_params, 'angles.itp')
write_dihedrals(atoms, dihedral_params, 'dihedrals.itp')

# STEP 8. Compose .itp files
compose_itp_files(['atomtypes.itp', 'moleculetype.itp', 'atoms.itp', 'bonds.itp', 'angles.itp', 'dihedrals.itp'], 'zif7.itp')

# STEP 9. Compose .itp files
with open('__tmp/atoms_zif7.pickle', 'wb') as handle:
    pickle.dump(atoms, handle, protocol=pickle.HIGHEST_PROTOCOL)







