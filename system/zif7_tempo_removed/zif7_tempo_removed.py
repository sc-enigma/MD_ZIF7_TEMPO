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

from zif7_params import get_zif7_params
from zif7_utils import transform_lattice

def find_loc(atoms, idx):
    r_best = np.zeros(3)
    dist_best = 0.0
    for theta in np.linspace(0, np.pi, 10):
        for phi in np.linspace(0, 2.0 * np.pi, 10):
            r = atoms[idx].r
            r[0] += 1.25 * np.sin(theta) * np.cos(phi)
            r[1] += 1.25 * np.sin(theta) * np.sin(phi)
            r[2] += 1.25 * np.cos(theta)
            
            dist = 1e10
            for atom_idx in range(len(atoms)):
                if atom_idx == idx:
                    continue
                delta = pow(sum(np.power(atoms[atom_idx].r - r, 2)), 0.5)
                dist = min(dist, delta)
            if dist > dist_best:
                dist_best = dist
                r_best = r
    return r_best

def add_OH(atoms, idx):
    l = len(atoms)
    atoms.append(Atom('OR', find_loc(atoms, idx)))
    atoms[idx].adjacency.append(l)
    atoms[l].adjacency.append(idx)
    atoms.append(Atom('HR', find_loc(atoms, l)))
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
    atoms.append(Atom('OR', find_loc(atoms, idx)))
    atoms[idx].adjacency.append(l)
    atoms[l].adjacency.append(idx)
    atoms.append(Atom('HR', find_loc(atoms, l)))
    atoms[l].adjacency.append(l + 1)
    atoms[l + 1].adjacency.append(l)
    atoms.append(Atom('HR', find_loc(atoms, l)))
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
    
#for
for removed_count in range(1, 4):
    # STEP 0. Load atoms   
    with open('../zif7_tempo/__tmp/atoms_tempo_np.pickle', 'rb') as handle:
        atoms_tempo_np = pickle.load(handle)
    with open('../zif7_tempo/__tmp/atoms_tempo_lp.pickle', 'rb') as handle:
        atoms_tempo_lp = pickle.load(handle)
    with open('../zif7/__tmp/atoms_zif7_np.pickle', 'rb') as handle:
        atoms_zif7_np = pickle.load(handle)
    with open('../zif7/__tmp/atoms_zif7_lp.pickle', 'rb') as handle:
        atoms_zif7_lp = pickle.load(handle)
    
    # STEP 1. Remove linkers    
    if removed_count == 1:
        ids_to_remove = [2275, 2293, 2285, 2262, 2263, 2264, 2268, 2269, 2290, 2291, 2286, 2287]
    if removed_count == 2:
        ids_to_remove = [2275, 2293, 2285, 2262, 2263, 2264, 2268, 2269, 2290, 2291, 2286, 2287,\
                         2208, 2190, 2200, 2177, 2202, 2201, 2179, 2178, 2205, 2206, 2183, 2184]
    if removed_count == 3:
        ids_to_remove = [2275, 2293, 2285, 2262, 2263, 2264, 2268, 2269, 2290, 2291, 2286, 2287,\
                         2208, 2190, 2200, 2177, 2202, 2201, 2179, 2178, 2205, 2206, 2183, 2184,\
                         2232, 2231, 2249, 2248, 2227, 2228, 2226, 2242, 2243, 2244, 2255, 2234]
    atoms_zif7_lp = remove_atoms(atoms_zif7_lp, ids_to_remove).copy()
    atoms_zif7_np = remove_atoms(atoms_zif7_np, ids_to_remove).copy()
    atoms_zif7_lp = remove_non_bonded_atoms(atoms_zif7_lp).copy()
    atoms_zif7_np = remove_non_bonded_atoms(atoms_zif7_np).copy()

    # STEP 2. Add OH and OH2
    if removed_count == 1:
        # 2270 2283
        atoms_zif7_np = add_OH(atoms_zif7_np,  2270)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2270)
        atoms_zif7_np = add_OH2(atoms_zif7_np, 2283)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2283)
    if removed_count == 2:
        # 2258 2271 2244 2224
        atoms_zif7_np = add_OH(atoms_zif7_np,  2258)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2258)
        atoms_zif7_np = add_OH2(atoms_zif7_np, 2271)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2271)
        atoms_zif7_np = add_OH(atoms_zif7_np,  2244)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2244)
        atoms_zif7_np = add_OH2(atoms_zif7_np, 2224)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2224)
    if removed_count == 3:
        # 2246 2259 2232 2218 2185 2198
        atoms_zif7_np = add_OH(atoms_zif7_np,  2246)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2246)
        atoms_zif7_np = add_OH2(atoms_zif7_np, 2259)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2259)
        atoms_zif7_np = add_OH(atoms_zif7_np,  2232)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2232)
        atoms_zif7_np = add_OH2(atoms_zif7_np, 2218)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2218)
        atoms_zif7_np = add_OH(atoms_zif7_np,  2185)
        atoms_zif7_lp = add_OH(atoms_zif7_lp,  2185)
        atoms_zif7_np = add_OH2(atoms_zif7_np, 2198)
        atoms_zif7_lp = add_OH2(atoms_zif7_lp, 2198)
        
    # STEP 3. Create .itp file
    mass, charge, bond_params, angle_params, dihedral_params = get_zif7_params()
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
    with open(f'removed_{removed_count}/__tmp/atoms_tempo_lp.pickle', 'wb') as handle:
        pickle.dump(atoms_tempo_lp, handle, protocol=pickle.HIGHEST_PROTOCOL)
    atoms_zif7_tempo_lp = atoms_zif7_lp.copy()
    atoms_zif7_tempo_lp.extend(atoms_tempo_lp.copy())
    write_gro_file(atoms_zif7_tempo_lp, f'removed_{removed_count}/zif7_tempo_lp_removed{removed_count}.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)

    # structure with nano pores
    a_np     = 24.60
    b_np     = 21.10
    c_np     = 15.70
    alpha_np = 84.80 / 180 * np.pi
    beta_np  = 86.10 / 180 * np.pi
    gamma_np = 128.00 / 180 * np.pi

    with open(f'removed_{removed_count}/__tmp/atoms_zif7_np.pickle', 'wb') as handle:
        pickle.dump(atoms_zif7_np, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(f'removed_{removed_count}/__tmp/atoms_tempo_np.pickle', 'wb') as handle:
        pickle.dump(atoms_tempo_np, handle, protocol=pickle.HIGHEST_PROTOCOL)
    atoms_zif7_tempo_np = atoms_zif7_np.copy()
    atoms_zif7_tempo_np.extend(atoms_tempo_np.copy())
    write_gro_file(atoms_zif7_tempo_np, f'removed_{removed_count}/zif7_tempo_np_removed{removed_count}.gro', a_np, b_np, c_np, alpha_np, beta_np, gamma_np, bounds_a, bounds_b, bounds_c)
