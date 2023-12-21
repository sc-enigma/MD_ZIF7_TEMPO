import sys
import numpy as np
import pickle
from numba import njit
from random import random

sys.path.append('../components')
sys.path.append('../zif7')
sys.path.append('../zif7_tempo')

from atom import Atom, mol2_to_atoms, remove_atoms, remove_non_bonded_atoms, copy_atom, count_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, compose_itp_files
from pbc import calculate_lattice_vectors

from zif7_params import get_zif7_params
from zif7_utils import transform_lattice

a     = 22.989
b     = 22.989
c     = 15.763
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 120.00 / 180 * np.pi
bounds_a, bounds_b, bounds_c = [1.0, 3.0], [1.0, 3.0], [1.0, 3.0]
    
@njit
def select_loc(atoms_r_internal, mol_r_internal, num_iter=100000):
    r_internal_best = np.zeros(3)
    dist_best = 0.0
    for iter in range(num_iter):
        r_internal = np.array([random() * 2.0, random() * 2.0, random() * 2.0])
        for i in range(len(atoms_r_internal)):
            dist = 1e10
            for j in range(len(mol_r_internal)):
                delta = np.abs(atoms_r_internal[i] - (mol_r_internal[j] + r_internal))
                delta[0] = min(delta[0], 2.0 - delta[0])
                delta[1] = min(delta[1], 2.0 - delta[1])
                delta[2] = min(delta[2], 2.0 - delta[2])
                dist = min(dist, pow(sum(np.power(delta, 2)), 0.5))
        if dist > dist_best:
            dist_best = dist
            r_internal_best = r_internal
        if dist_best > 1.4:
            break
    return r_internal_best

alpha_rev = np.arccos((np.cos(beta) * np.cos(gamma) - np.cos(alpha)) / np.sin(beta) / np.sin(gamma))
transformation = np.array([\
    [a,   b*np.cos(gamma), c*np.cos(beta)], \
    [0.0, b*np.sin(gamma), -1.0*c*np.sin(beta)*np.cos(alpha_rev)], \
    [0.0, 0.0,             c*np.sin(beta)*np.sin(alpha_rev)
]])
vec_a, vec_b, vec_c = calculate_lattice_vectors(a, b, c, alpha, beta, gamma)
transformation_inv = np.linalg.inv(transformation)
    
    
for removed_count in range(4):
    if removed_count == 0:
        with open('../zif7_tempo/__tmp/atoms_tempo_lp.pickle', 'rb') as handle:
            atoms_tempo = pickle.load(handle)
        with open('../zif7/__tmp/atoms_zif7_lp.pickle', 'rb') as handle:
            atoms_zif7 = pickle.load(handle)
    else:
        with open(f'../zif7_tempo_removed/removed_{removed_count}/__tmp/atoms_tempo_lp.pickle', 'rb') as handle:
            atoms_tempo = pickle.load(handle)
        with open(f'../zif7_tempo_removed/removed_{removed_count}/__tmp/atoms_zif7_lp.pickle', 'rb') as handle:
            atoms_zif7 = pickle.load(handle)
    atoms_zif7_tempo = atoms_zif7.copy()
    atoms_zif7_tempo.extend(atoms_tempo.copy())
    
    mol_count = 50
    for mol_idx in range(mol_count):
        print(mol_idx)
        for atom_idx in range(len(atoms_zif7_tempo)):
            atoms_zif7_tempo[atom_idx].r_internal = transformation_inv.dot(np.transpose(atoms_zif7_tempo[atom_idx].r))
                                                                                            
        atoms_r_internal = np.array([atom.r_internal for atom in atoms_zif7_tempo])
        mol_r_internal = [transformation_inv.dot(np.transpose(r)) for r in [\
            np.array([1.160,  0.000, 0.000]),
            np.array([0.000,  0.000, 0.000]),
            np.array([-1.160, 0.000, 0.000]),
            np.array([-0.990, 0.000, 0.000]),
            np.array([0.990,  0.000, 0.000]),
            ]]
        mol_r_internal = np.array(mol_r_internal)
        shift = select_loc(atoms_r_internal, mol_r_internal)

        mol_r = np.empty((len(mol_r_internal), 3))
        for j in range(len(mol_r_internal)):
            mol_r_internal[j] += shift
            mol_r[j] = mol_r_internal[j][0] * vec_a + mol_r_internal[j][1] * vec_b + mol_r_internal[j][2] * vec_c
        
        mol_atoms = [\
            Atom('O1', mol_r[0]),
            Atom('C',  mol_r[1]),
            Atom('O2', mol_r[2]),
            Atom('M1', mol_r[3]),
            Atom('M2', mol_r[4])]
            
        mol_atoms[0].resid = 'CO2'
        mol_atoms[1].resid = 'CO2'
        mol_atoms[2].resid = 'CO2'
        mol_atoms[3].resid = 'CO2'
        mol_atoms[4].resid = 'CO2'
        
        mol_atoms[0].resid_idx = mol_idx + 1
        mol_atoms[1].resid_idx = mol_idx + 1
        mol_atoms[2].resid_idx = mol_idx + 1
        mol_atoms[3].resid_idx = mol_idx + 1
        mol_atoms[4].resid_idx = mol_idx + 1
        
        mol_atoms[0].atom_idx = len(atoms_zif7_tempo)
        mol_atoms[1].atom_idx = len(atoms_zif7_tempo) + 1
        mol_atoms[2].atom_idx = len(atoms_zif7_tempo) + 2
        mol_atoms[3].atom_idx = len(atoms_zif7_tempo) + 3
        mol_atoms[4].atom_idx = len(atoms_zif7_tempo) + 4
         
        atoms_zif7_tempo.append(copy_atom(mol_atoms[0]))
        atoms_zif7_tempo.append(copy_atom(mol_atoms[1]))
        atoms_zif7_tempo.append(copy_atom(mol_atoms[2]))
        atoms_zif7_tempo.append(copy_atom(mol_atoms[3]))
        atoms_zif7_tempo.append(copy_atom(mol_atoms[4]))
    
    write_gro_file(atoms_zif7_tempo, f'zif7_tempo_CO2_lp_removed{removed_count}.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
