import numpy as np

from atom import Atom

def getAllBonds(atoms):
    bonds = []
    for atom_idx in range(len(atoms)):
        adjacency = atoms[atom_idx].adjacency
        for adj in adjacency:
            bond = np.sort([atom_idx, adj])
            
            counts = np.unique(bond, return_counts=True)[1]
            correct = len(np.unique(counts)) == 1

            if correct:
                bonds.append(bond)
    bonds = np.array(bonds)
    return np.unique(bonds, axis=0)

def getAllAngles(atoms):
    angles = []
    for atom_idx in range(len(atoms)):
        adjacency1 = atoms[atom_idx].adjacency
        for adj1 in adjacency1:
            adjacency2 = atoms[adj1].adjacency
            for adj2 in adjacency2:
                angle = np.array([atom_idx, adj1, adj2])
                if angle[0] > angle[2]:
                    angle = np.flip(angle)
            
                counts = np.unique(angle, return_counts=True)[1]
                correct = len(np.unique(counts)) == 1

                if correct:
                    angles.append(angle)
    
    angles = np.array(angles)
    return np.unique(angles, axis=0)

def getAllDihedrals(atoms):
    dihedrals = []
    for atom_idx in range(len(atoms)):
        adjacency1 = atoms[atom_idx].adjacency
        for adj1 in adjacency1:
            adjacency2 = atoms[adj1].adjacency
            for adj2 in adjacency2:
                adjacency3 = atoms[adj2].adjacency
                for adj3 in adjacency3:
                    dihedral = np.array([atom_idx, adj1, adj2, adj3])
                    if dihedral[0] > dihedral[3]:
                        dihedral = np.flip(dihedral)
                    
                    counts = np.unique(dihedral, return_counts=True)[1]
                    correct = len(np.unique(counts)) == 1
                    correct &= np.unique(counts)[0] == 1

                    if correct:
                        dihedrals.append(dihedral)
                        
    dihedrals = np.array(dihedrals)
    return np.unique(dihedrals, axis=0)
