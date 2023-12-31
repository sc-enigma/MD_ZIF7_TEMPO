import numpy as np
import matplotlib.pyplot as plt

class Atom:
    def __init__(self, name, r):
        self.name = name        
        self.r = np.array(r)
        
        # topology properties
        self.resid_idx = 0
        self.resid = 'UNK'
        self.atom_idx = 0
        self.atom_type = 'unk_type'
        
        # internal properties
        self.r_internal = np.zeros(3)
        self.adjacency = []

def copy_atom(atom):
    a = Atom(atom.name, np.array(atom.r))
    a.resid_idx = atom.resid_idx
    a.resid = atom.resid
    a.atom_idx = atom.atom_idx
    a.atom_type = atom.atom_type
    a.r_internal = np.array(a.r_internal)
    a.adjacency = atom.adjacency.copy()
    return a
        
def mol2_to_atoms(data):
    atoms = []
    
    data_atoms = data['ATOM']
    for atom_idx in range(len(data_atoms)):
        split = data_atoms[atom_idx].split()
        if len(split) < 5:
            continue
        name = split[1]
        r = [float(c) for c in split[2:5]]
        atoms.append(Atom(name, r))
    
    if 'BOND' in data.keys():
        data_bonds = data['BOND']
        for bond_idx in range(len(data_bonds)):
            split = data_bonds[bond_idx].split()
            connected = [int(split[1]) - 1, int(split[2]) - 1]
            atoms[connected[0]].adjacency.append(connected[1])
            atoms[connected[1]].adjacency.append(connected[0])
        
    return atoms
    
def count_atoms(atoms):
    count = {}
    for atom_idx in range(len(atoms)):
        atom_type = atoms[atom_idx].atom_type
        if atom_type in count.keys():
            count[atom_type] += 1
        else:
            count[atom_type] = 1
    for atom_type in count.keys():
        print(f'{atom_type}    {count[atom_type]}')
        
def dump_atoms(atoms, dim_1=0, dim_2=1, dumpBonds=False):
    points_single = []
    point_pairs = []
    for atom_idx in range(len(atoms)):
        points_single.append(atoms[atom_idx].r)
        for adj_idx in atoms[atom_idx].adjacency:
            if atoms[atom_idx].alive and atoms[adj_idx].alive:
                point_pairs.append([atoms[atom_idx].r, atoms[adj_idx].r])
        
    points_single = np.transpose(np.array(points_single))
    s = np.ones(len(points_single[0])) * 0.5
    plt.scatter(points_single[dim_1], points_single[dim_2], s=s)
    
    if dumpBonds:
        for pair_idx in range(len(point_pairs)):
            pair = point_pairs[pair_idx]
            plt.plot([pair[0][dim_1], pair[1][dim_1]], [pair[0][dim_2], pair[1][dim_2]])        
    plt.show()

def remove_atoms(atoms, ids_to_remove):
    if len(ids_to_remove) == 0:
        return atoms
    shift_ids = {}
    removed_count = 0
    for atom_idx in range(len(atoms)):
        if atom_idx in ids_to_remove:
            shift_ids[atom_idx] = -1
            removed_count += 1
        else:
            shift_ids[atom_idx] = atom_idx - removed_count
            
    atoms = [copy_atom(atoms[atom_idx]) for atom_idx in range(len(atoms)) if not (atom_idx in ids_to_remove)]

    for atom_idx in range(len(atoms)):
        adjacency = atoms[atom_idx].adjacency
        for adj_idx in range(len(adjacency)):
            atoms[atom_idx].adjacency[adj_idx] = shift_ids[adjacency[adj_idx]]
            
        while -1 in atoms[atom_idx].adjacency:
            atoms[atom_idx].adjacency.remove(-1)
    return atoms
    
def remove_non_bonded_atoms(atoms):
    ids_to_remove = []
    for atom_idx in range(len(atoms)):
        if len(atoms[atom_idx].adjacency) == 0:
            ids_to_remove.append(atom_idx)
    return remove_atoms(atoms, ids_to_remove)