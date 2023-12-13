import numpy as np
import pickle

from atom import Atom

# Calculate crystal lattice vectors
def calculate_lattice_vectors(a, b, c,
                              alpha, beta, gamma):
    alpha_rev = np.arccos((np.cos(beta) * np.cos(gamma) - np.cos(alpha)) / np.sin(beta) / np.sin(gamma))

    transformation = np.array([\
        [a,     b*np.cos(gamma),   c*np.cos(beta)                       ], \
        [0.0,   b*np.sin(gamma),   -1.0*c*np.sin(beta)*np.cos(alpha_rev)], \
        [0.0,   0.0,               c*np.sin(beta)*np.sin(alpha_rev)     ]
        ])

    vec_a = transformation.dot(np.transpose([1.0, 0.0, 0.0]))
    vec_b = transformation.dot(np.transpose([0.0, 1.0, 0.0]))
    vec_c = transformation.dot(np.transpose([0.0, 0.0, 1.0]))
    
    # print(f'VECTOR A = {vec_a}')
    # print(f'VECTOR B = {vec_b}')
    # print(f'VECTOR C = {vec_c}')
    return vec_a, vec_b, vec_c

def compare_periodic(r1, r2, p_a, p_b, p_c):
    delta = r1 - r2 + np.array([p_a, p_b, p_c]) * 10.0
    tol = 1e-3

    delta[0] = delta[0] % p_a
    delta[1] = delta[1] % p_b
    delta[2] = delta[2] % p_c
    
    delta[0] = min(delta[0], p_a - delta[0])
    delta[1] = min(delta[1], p_b - delta[1])
    delta[2] = min(delta[2], p_c - delta[2])
    return abs(delta[0]) < tol and abs(delta[1]) < tol and abs(delta[2]) < tol

def apply_PBC(atoms,
              a, b, c,
              alpha, beta, gamma,
              bounds_a, bounds_b, bounds_c):
    # Calculate internal coordinates
    alpha_rev = np.arccos((np.cos(beta) * np.cos(gamma) - np.cos(alpha)) / np.sin(beta) / np.sin(gamma))

    transformation = np.array([\
        [a,   b*np.cos(gamma), c*np.cos(beta)], \
        [0.0, b*np.sin(gamma), -1.0*c*np.sin(beta)*np.cos(alpha_rev)], \
        [0.0, 0.0,             c*np.sin(beta)*np.sin(alpha_rev)
    ]])
    
    transformation_inv = np.linalg.inv(transformation)
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].r_internal = transformation_inv.dot(np.transpose(atoms[atom_idx].r))

    # to create unique.npy and adjacency_data.pickle decomment code below
    '''tol_internal = 1e-4
    unique = np.empty(len(atoms))
    # define what atoms are unique and what are copies
    for atom_idx in range(len(atoms)):
        r_internal = atoms[atom_idx].r_internal
        unique[atom_idx] = (\
            r_internal[0] > bounds_a[0] + tol_internal and r_internal[0] < bounds_a[1] + tol_internal and
            r_internal[1] > bounds_b[0] + tol_internal and r_internal[1] < bounds_b[1] + tol_internal and
            r_internal[2] > bounds_c[0] + tol_internal and r_internal[2] < bounds_c[1] + tol_internal)
    
    # remove duplicates within unique atoms
    for atom_idx1 in range(len(atoms)):
        if not unique[atom_idx1]:
            continue
        for atom_idx2 in range(atom_idx1 + 1, len(atoms)):
            if not unique[atom_idx2]:
                continue
            if sum(np.abs(atoms[atom_idx1].r - atoms[atom_idx2].r)) < 1e-3:
                unique[atom_idx2] = False
                for adj in atoms[atom_idx2].adjacency:
                    atoms[atom_idx1].adjacency.append(adj)
                    atoms[atom_idx1].adjacency = np.unique(atoms[atom_idx1].adjacency).tolist()

    # correct adjacency: remove ids of copies by ids of their originals
    tol_internal = 1e-4
    p_a, p_b, p_c = bounds_a[1] - bounds_a[0], bounds_b[1] - bounds_b[0], bounds_c[1] - bounds_c[0]
    for atom_idx in range(len(atoms)):
        if not unique[atom_idx]:
            continue
        adjacency = atoms[atom_idx].adjacency
        for adj_idx in range(len(adjacency)):
            if unique[adjacency[adj_idx]]:
                continue
            r_internal = atoms[adjacency[adj_idx]].r_internal
            
            orig_cnt = 0
            orig_idx = 0
            for cmp_idx in range(len(atoms)):
                if not unique[cmp_idx]:
                    continue
                cmp_r_internal = atoms[cmp_idx].r_internal
                if compare_periodic(r_internal, cmp_r_internal, p_a, p_b, p_c):
                    orig_idx = cmp_idx
                    orig_cnt += 1
            
            atoms[atom_idx].adjacency[adj_idx] = orig_idx
            if orig_cnt != 1:
                print('ERROR', orig_cnt, r_internal)
    
    np.save('__tmp/unique.npy', unique)
    adjacency_data = [atom.adjacency for atom in atoms]
    with open('__tmp/adjacency_data.pickle', 'wb') as handle:
        pickle.dump(adjacency_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    '''
    unique = np.load('__tmp/unique.npy')
    with open('__tmp/adjacency_data.pickle', 'rb') as handle:
        adjacency_data = pickle.load(handle)
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].adjacency = np.unique(adjacency_data[atom_idx].copy()).tolist()

    # remove copies from atoms
    shift_ids = {}
    alive_count = 0
    for atom_idx in range(len(atoms)):
        if unique[atom_idx]:
            shift_ids[atom_idx] = alive_count
            alive_count += 1
    atoms = [atoms[atom_idx] for atom_idx in range(len(atoms)) if unique[atom_idx] > 0.0]
    
    # update adjacency after copies removal
    for atom_idx in range(len(atoms)):
        adjacency = atoms[atom_idx].adjacency
        for adj_idx in range(len(adjacency)):
            atoms[atom_idx].adjacency[adj_idx] = shift_ids[adjacency[adj_idx]]
            if adjacency[adj_idx] == atom_idx or adjacency[adj_idx] >= len(atoms):
                print('ERROR')
    
    # shift atoms
    vec_a, vec_b, vec_c = calculate_lattice_vectors(a, b, c, alpha, beta, gamma)
    shift = bounds_a[0] * vec_a + bounds_b[0] * vec_b + bounds_c[0] * vec_c
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].r -= shift
    return atoms
