import sys
import numpy as np
import pickle

sys.path.append('../components')
sys.path.append('../zif7')
from atom import Atom, mol2_to_atoms
from read_utils import read_mol2_file
from write_utils import write_gro_file, write_mol2_file
from pbc import calculate_lattice_vectors

from zif7_utils import transform_lattice

with open('../zif7/__tmp/atoms_zif7_np.pickle', 'rb') as handle:
    atoms_zif7_np = pickle.load(handle)
with open('../zif7/__tmp/atoms_zif7_lp.pickle', 'rb') as handle:
    atoms_zif7_lp = pickle.load(handle)
atoms_tempo = mol2_to_atoms(read_mol2_file('../tempo/tempo.mol2'))
for atom_idx in range(len(atoms_tempo)):
    atoms_tempo[atom_idx].resid_idx = 1
    atoms_tempo[atom_idx].resid = 'TMP'
    atoms_tempo[atom_idx].atom_idx = atom_idx
 
def put_tempo_in_lattice(atoms_zif7, atoms_tempo, shift=np.zeros(3)):
    # Calculate cell center
    cell_center = np.zeros(3)
    for atom_idx in [2275, 2294, 2255, 2235, 2190, 2209]:
        cell_center += atoms_zif7[atom_idx].r / 6.0

    # Calculate TEMPO center
    tempo_center = np.zeros(3)
    for atom in atoms_tempo:
        tempo_center += atom.r / len(atoms_tempo)
        
    # Translate TEMPO atoms inside cell
    for atom_idx in range(len(atoms_tempo)):
        atoms_tempo[atom_idx].r += cell_center - tempo_center + shift

    # Calculate minimal distance from TEMPO to cell
    min_dist = 1e10
    for atom_tempo in atoms_tempo:
        for atom_zif7 in atoms_zif7:
            dist = pow(sum(np.power(atom_tempo.r - atom_zif7.r, 2)), 0.5)
            min_dist = min(min_dist, dist)
    print(min_dist)
    return atoms_zif7.copy(), atoms_tempo.copy()
    
# structure with large pores
a     = 22.989
b     = 22.989
c     = 15.763
alpha = 90.00 / 180 * np.pi
beta  = 90.00 / 180 * np.pi
gamma = 120.00 / 180 * np.pi
bounds_a, bounds_b, bounds_c = [1.0, 3.0], [1.0, 3.0], [1.0, 3.0]

# Save .gro file
atoms_zif7_lp, atoms_tempo = put_tempo_in_lattice(atoms_zif7_lp, atoms_tempo)
atoms_zif7_tempo_lp = atoms_zif7_lp.copy()
atoms_zif7_tempo_lp.extend(atoms_tempo.copy())
write_gro_file(atoms_zif7_tempo_lp, 'zif7_tempo_lp.gro', a, b, c, alpha, beta, gamma, bounds_a, bounds_b, bounds_c)

# structure with nano pores
a_np     = 24.60
b_np     = 21.10
c_np     = 15.70
alpha_np = 84.80 / 180 * np.pi
beta_np  = 86.10 / 180 * np.pi
gamma_np = 128.00 / 180 * np.pi

# Save .gro file
atoms_zif7_np, atoms_tempo = put_tempo_in_lattice(atoms_zif7_np, atoms_tempo, np.array([0.7, 0.2, 0.1]))
atoms_zif7_tempo_np = atoms_zif7_np.copy()
atoms_zif7_tempo_np.extend(atoms_tempo.copy())
write_gro_file(atoms_zif7_tempo_np, 'zif7_tempo_np.gro', a_np, b_np, c_np, alpha_np, beta_np, gamma_np, bounds_a, bounds_b, bounds_c)
