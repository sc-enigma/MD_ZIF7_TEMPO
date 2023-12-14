import numpy as np

from atom import Atom
from pbc import calculate_lattice_vectors

def remove_oxygen(atoms):
    shift_ids = {}
    oxygen_count = 0
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'O':
            shift_ids[atom_idx] = -1
            oxygen_count += 1
        else:
            shift_ids[atom_idx] = atom_idx - oxygen_count
    atoms = [atom for atom in atoms if atom.name[0] != 'O']

    for atom_idx in range(len(atoms)):
        adjacency = atoms[atom_idx].adjacency
        for adj_idx in range(len(adjacency)):
            atoms[atom_idx].adjacency[adj_idx] = shift_ids[adjacency[adj_idx]]
            
        if -1 in atoms[atom_idx].adjacency:
            atoms[atom_idx].adjacency.remove(-1)
    return atoms
    
def define_zif7_atom_types(atoms):
    def getNeighbours(atom_idx):
        return [atoms[adj_idx] for adj_idx in atoms[atom_idx].adjacency]
    
    def getNeighbourElems(atom_idx):
        return np.sort([atom.name[0] for atom in getNeighbours(atom_idx)])
    
    def getElem(atom_idx):
        return atoms[atom_idx].name[0]
        
    # define N: intraFF_N0
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'N':
            atoms[atom_idx].atom_type = 'intraFF_N'
            
    # define Zn: intraFF_Zn
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].name[0] == 'Z':
            atoms[atom_idx].atom_type = 'intraFF_Zn'
  
    # define C1: intraFF_C1
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'C' and atoms[atom_idx].atom_type == '':
            if (np.array_equal(getNeighbourElems(atom_idx), ['H', 'N', 'N'])):
                atoms[atom_idx].atom_type = 'intraFF_C1'
    
    # define C2: intraFF_C2
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'C' and atoms[atom_idx].atom_type == '':
            if (np.array_equal(getNeighbourElems(atom_idx), ['C', 'C', 'N'])):
                atoms[atom_idx].atom_type = 'intraFF_C2'
                
    # define C5: intraFF_C5
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'C' and atoms[atom_idx].atom_type == '':
            for neighbour in getNeighbours(atom_idx):
                if neighbour.atom_type == 'intraFF_C2':
                    atoms[atom_idx].atom_type = 'intraFF_C5'
                    break
                
    # define C6: intraFF_C6
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'C' and atoms[atom_idx].atom_type == '':
            for neighbour in getNeighbours(atom_idx):
                if neighbour.atom_type == 'intraFF_C5':
                    atoms[atom_idx].atom_type = 'intraFF_C6'
                    break
                
    # define H1, H5, H6: intraFF_H1, intraFF_H5, intraFF_H6
    for atom_idx in range(len(atoms)):
        if getElem(atom_idx) == 'H':
            neighbour = getNeighbours(atom_idx)[0]
            if neighbour.atom_type == 'intraFF_C1':
                atoms[atom_idx].atom_type = 'intraFF_H1'
            if neighbour.atom_type == 'intraFF_C5':
                atoms[atom_idx].atom_type = 'intraFF_H5'
            if neighbour.atom_type == 'intraFF_C6':
                atoms[atom_idx].atom_type = 'intraFF_H6'
            
    # check that all atom types defined
    for atom_idx in range(len(atoms)):
        if atoms[atom_idx].atom_type == '':
            print('ERROR')
            
            print(getElem(atom_idx), getNeighbourElems(atom_idx))
            
    return atoms

def define_zif7_atom_names(atoms):  
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].name = atoms[atom_idx].atom_type.replace('intraFF_', '')
        atoms[atom_idx].resid_idx = 1
        atoms[atom_idx].resid = 'ZIF'
        atoms[atom_idx].atom_idx = atom_idx
        
    return atoms

def transform_lattice(atoms,
                      a, b, c, alpha, beta, gamma,
                      a_new, b_new, c_new, alpha_new, beta_new, gamma_new):
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
    
    # Calculate internal coordinates
    vec_a_new, vec_b_new, vec_c_new = calculate_lattice_vectors(a_new, b_new, c_new, alpha_new, beta_new, gamma_new)
    for atom_idx in range(len(atoms)):
        atoms[atom_idx].r = \
                          atoms[atom_idx].r_internal[0] * vec_a_new +\
                          atoms[atom_idx].r_internal[1] * vec_b_new +\
                          atoms[atom_idx].r_internal[2] * vec_c_new
    return atoms



