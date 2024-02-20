import sys
import numpy as np
import pickle

sys.path.append('../components')
sys.path.append('../zif7')
sys.path.append('../zif7_tempo')
from atom import Atom, mol2_to_atoms, remove_atoms, remove_non_bonded_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from write_utils import write_atoms, write_bonds, write_angles, write_dihedrals, compose_itp_files
from pbc import calculate_lattice_vectors

from zif7_params_lp import get_zif7_params_lp
from zif7_utils import transform_lattice

def find_loc(atoms, idx):
    r_best = np.zeros(3)
    dist_best = 0.0
    for theta in np.linspace(0, np.pi, 10):
        for phi in np.linspace(0, 2.0 * np.pi, 10):
            r = np.copy(atoms[idx].r)
            r[0] += 1.25 * np.sin(theta) * np.cos(phi)
            r[1] += 1.25 * np.sin(theta) * np.sin(phi)
            r[2] += 1.25 * np.cos(theta)
            
            dist = 1e10
            for atom_idx in range(len(atoms)):
                if atom_idx == idx:
                    continue
                delta = pow(sum(np.power(np.copy(atoms[atom_idx].r) - np.copy(r), 2)), 0.5)
                dist = min(dist, delta)
            if dist > dist_best:
                dist_best = dist
                r_best = np.copy(r)
    return r_best

def add_OH(atoms, idx):
    l = len(atoms)
    atoms.append(Atom('OR', np.copy(find_loc(atoms, idx))))
    atoms[idx].adjacency.append(l)
    atoms[l].adjacency.append(idx)
    atoms.append(Atom('HR', np.copy(find_loc(atoms, l))))
    atoms[l].adjacency.append(l + 1)
    atoms[l + 1].adjacency.append(l)
    
    atoms[l].resid_idx = 1
    atoms[l + 1].resid_idx = 1
    atoms[l].resid = 'ZIF'
    atoms[l + 1].resid = 'ZIF'
    atoms[l].atom_idx = 1
    atoms[l + 1].atom_idx = 1
    atoms[l].atom_type = 'repl_O'
    atoms[l + 1].atom_type = 'repl_H'
    
    return atoms
    
def add_OH2(atoms, idx):
    l = len(atoms)
    atoms.append(Atom('OR', np.copy(find_loc(atoms, idx))))
    atoms[idx].adjacency.append(l)
    atoms[l].adjacency.append(idx)
    atoms.append(Atom('HR', np.copy(find_loc(atoms, l))))
    atoms[l].adjacency.append(l + 1)
    atoms[l + 1].adjacency.append(l)
    atoms.append(Atom('HR', np.copy(find_loc(atoms, l))))
    atoms[l].adjacency.append(l + 2)
    atoms[l + 2].adjacency.append(l)
    
    atoms[l].resid_idx = 1
    atoms[l + 1].resid_idx = 1
    atoms[l + 2].resid_idx = 1
    atoms[l].resid = 'ZIF'
    atoms[l + 1].resid = 'ZIF'
    atoms[l + 2].resid = 'ZIF'
    atoms[l].atom_idx = 1
    atoms[l + 1].atom_idx = 1
    atoms[l + 2].atom_idx = 1
    atoms[l].atom_type = 'repl_O'
    atoms[l + 1].atom_type = 'repl_H'
    atoms[l + 2].atom_type = 'repl_H'
    
    return atoms
    
for removed_count in range(1, 4):
    # STEP 0. Load atoms   
    with open('../zif7_tempo/__tmp/atoms_tempo_lpA.pickle', 'rb') as handle:
        atoms_tempo_lp = pickle.load(handle)
    with open('../zif7/__tmp/atoms_zif7_lp.pickle', 'rb') as handle:
        atoms_zif7_lp = pickle.load(handle)
    
    # STEP 1. Remove linkers
    # z < 17 and z > 9 and y < 20 and resname ZIF
    if removed_count == 1:
        ids_to_remove = [2350, 2364, 2389, 2351, 2381, 2383, 2382, 2387, 2386, 2357, 2358, 2352, 2353]
    if removed_count == 2:
        ids_to_remove = [2350, 2364, 2389, 2351, 2381, 2383, 2382, 2387, 2386, 2357, 2358, 2352, 2353,\
                         2454, 2468, 2493, 2455, 2485, 2457, 2456, 2486, 2487, 2461, 2462, 2490, 2491]
    if removed_count == 3:
        ids_to_remove = [2350, 2364, 2389, 2351, 2381, 2383, 2382, 2387, 2386, 2357, 2358, 2352, 2353,\
                         2454, 2468, 2493, 2455, 2485, 2457, 2456, 2486, 2487, 2461, 2462, 2490, 2491,\
                         2415, 2428, 2442, 2407, 2408, 2409, 2412, 2413, 2435, 2436, 2430, 2429, 2431]
    atoms_zif7_lp = remove_atoms(atoms_zif7_lp, ids_to_remove)
    
    # STEP 2. Add OH and OH2
    if removed_count == 1:
        pass
        # 2270 2283
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2378)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2358)
    if removed_count == 2:
        # 2258 2271 2244 2224
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2469)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2449)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2378)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2358)
    if removed_count == 3:
        # 2246 2259 2232 2218 2185 2198
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2456)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2436)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2417)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2398)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2378)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2358)
        
    atoms_zif7_lp = remove_non_bonded_atoms(atoms_zif7_lp)
        
    # STEP 3. Create .itp file
    mass, charge, bond_params, angle_params, dihedral_params = get_zif7_params_lp()
    write_atoms(atoms_zif7_lp, charge, mass, 'ZIF', f'removed_{removed_count}/atoms.itp')
    write_bonds(atoms_zif7_lp, bond_params,         f'removed_{removed_count}/bonds.itp')
    write_angles(atoms_zif7_lp, angle_params,       f'removed_{removed_count}/angles.itp')
    write_dihedrals(atoms_zif7_lp, dihedral_params, f'removed_{removed_count}/dihedrals.itp')
    compose_itp_files(['atomtypes.itp',
                       'moleculetype.itp',
                       f'removed_{removed_count}/atoms.itp',
                       f'removed_{removed_count}/bonds.itp',
                       f'removed_{removed_count}/angles.itp',
                       f'removed_{removed_count}/dihedrals.itp'],
                       filename=f'removed_{removed_count}/zif7.itp')

    # STEP 4. Save .gro files
    # structure with large pores
    a     = 22.989
    b     = 22.989
    c     = 15.763
    alpha = 90.00 / 180 * np.pi
    beta  = 90.00 / 180 * np.pi
    gamma = 120.00 / 180 * np.pi
    bounds_a, bounds_b, bounds_c = [1.0, 3.0], [1.0, 3.0], [1.0, 3.0]

    with open(f'removed_{removed_count}/__tmp/atoms_zif7_lp.pickle', 'wb') as handle:
        pickle.dump(atoms_zif7_lp, handle, protocol=pickle.HIGHEST_PROTOCOL)
    atoms_zif7_tempo_lp = atoms_zif7_lp.copy()
    atoms_zif7_tempo_lp.extend(atoms_tempo_lp.copy())
    write_gro_file(atoms_zif7_tempo_lp, f'removed_{removed_count}/zif7_tempo_lp_removed{removed_count}.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)
    
    print(removed_count)