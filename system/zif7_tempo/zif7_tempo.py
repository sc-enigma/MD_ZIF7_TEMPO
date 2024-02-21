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

atoms_tempo = mol2_to_atoms(read_mol2_file('../tempo/tempo.mol2'))
for atom_idx in range(len(atoms_tempo)):
    atoms_tempo[atom_idx].resid_idx = 1
    atoms_tempo[atom_idx].resid = 'TMP'
    atoms_tempo[atom_idx].atom_idx = atom_idx
 
def put_tempo_in_lattice(atoms_zif7, atoms_tempo, pore_idx):
    def calc_dist(atoms_zif7, atoms_tempo, shift):
        min_dist = 1e10
        for atom_tempo in atoms_tempo:
            for atom_zif7 in atoms_zif7:
                dist = pow(sum(np.power(atom_tempo.r + shift - atom_zif7.r, 2)), 0.5)
                min_dist = min(min_dist, dist)
        return min_dist
    # Calculate cell center
    cell_center = np.zeros(3)
    for atom_idx in pore_idx:
        cell_center += atoms_zif7[atom_idx].r / len(pore_idx)

    # Calculate TEMPO center
    tempo_center = np.zeros(3)
    for atom in atoms_tempo:
        tempo_center += atom.r / len(atoms_tempo)
        
    # Translate TEMPO atoms inside cell
    for atom_idx in range(len(atoms_tempo)):
        atoms_tempo[atom_idx].r += cell_center - tempo_center

    best_dist = 0.0
    best_shift = np.zeros(3)
    for dx in np.linspace(-2.0, 2.0, 5):
        for dy in np.linspace(-2.0, 2.0, 5):
            for dz in np.linspace(-2.0, 2.0, 5):
                shift = np.array([dx, dy, dz])
                dist = calc_dist(atoms_zif7, atoms_tempo, shift)
                print(dist)
                if dist > best_dist:
                    best_dist = dist
                    best_shift = shift
    for atom_idx in range(len(atoms_tempo)):
        atoms_tempo[atom_idx].r += best_shift

    # Calculate minimal distance from TEMPO to cell
    print(calc_dist(atoms_zif7, atoms_tempo, np.zeros(3)), best_shift)
    return atoms_zif7.copy(), atoms_tempo.copy()

# STEP 1. Structure with narrow pores. Pore A ________________________________________________________________________
a_np     = 22.989
b_np     = 22.989
c_np     = 15.763
alpha_np = 90.00 / 180 * np.pi
beta_np  = 90.00 / 180 * np.pi
gamma_np = 120.00 / 180 * np.pi
bounds_a, bounds_b, bounds_c = [1.0, 3.0], [1.0, 3.0], [1.0, 3.0]

# Load data
with open('../zif7/__tmp/atoms_zif7_np.pickle', 'rb') as handle:
    atoms_zif7_np = pickle.load(handle)
    
# Save .gro file
atoms_zif7_np, atoms_tempo = put_tempo_in_lattice(atoms_zif7_np, atoms_tempo, [2362, 2468])
with open('__tmp/atoms_zif7_npA.pickle', 'wb') as handle:
    pickle.dump(atoms_zif7_np, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('__tmp/atoms_tempo_npA.pickle', 'wb') as handle:
    pickle.dump(atoms_tempo, handle, protocol=pickle.HIGHEST_PROTOCOL)
atoms_zif7_tempo_npA = atoms_zif7_np.copy()
atoms_zif7_tempo_npA.extend(atoms_tempo.copy())
write_gro_file(atoms_zif7_tempo_npA, 'zif7_tempo_npA.gro', a_np, b_np, c_np, alpha_np, beta_np, gamma_np, bounds_a, bounds_b, bounds_c)

# STEP 2. Structure with narrow pores. Pore B ________________________________________________________________________
# Load data
with open('../zif7/__tmp/atoms_zif7_np.pickle', 'rb') as handle:
    atoms_zif7_np = pickle.load(handle)
    
# Save .gro file
atoms_zif7_np, atoms_tempo = put_tempo_in_lattice(atoms_zif7_np, atoms_tempo, [1783, 2465])
with open('__tmp/atoms_zif7_npB.pickle', 'wb') as handle:
    pickle.dump(atoms_zif7_np, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('__tmp/atoms_tempo_npB.pickle', 'wb') as handle:
    pickle.dump(atoms_tempo, handle, protocol=pickle.HIGHEST_PROTOCOL)
atoms_zif7_tempo_npB = atoms_zif7_np.copy()
atoms_zif7_tempo_npB.extend(atoms_tempo.copy())
write_gro_file(atoms_zif7_tempo_npB, 'zif7_tempo_npB.gro', a_np, b_np, c_np, alpha_np, beta_np, gamma_np, bounds_a, bounds_b, bounds_c)

# STEP 3. Structure with large pores. Pore A ________________________________________________________________________
a_lp     = 23.948
b_lp     = 21.354
c_lp     = 16.349
alpha_lp = 90.28 / 180 * np.pi
beta_lp  = 93.28 / 180 * np.pi
gamma_lp = 108.41 / 180 * np.pi

# Load data
with open('../zif7/__tmp/atoms_zif7_lp.pickle', 'rb') as handle:
    atoms_zif7_lp = pickle.load(handle)

# Save .gro file
atoms_zif7_lp, atoms_tempo = put_tempo_in_lattice(atoms_zif7_lp, atoms_tempo, [2362, 2468])
with open('__tmp/atoms_zif7_lpA.pickle', 'wb') as handle:
    pickle.dump(atoms_zif7_lp, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('__tmp/atoms_tempo_lpA.pickle', 'wb') as handle:
    pickle.dump(atoms_tempo, handle, protocol=pickle.HIGHEST_PROTOCOL)
atoms_zif7_tempo_lpA = atoms_zif7_lp.copy()
atoms_zif7_tempo_lpA.extend(atoms_tempo.copy())
write_gro_file(atoms_zif7_tempo_lpA, 'zif7_tempo_lpA.gro', a_lp, b_lp, c_lp, alpha_lp, beta_lp, gamma_lp, bounds_a, bounds_b, bounds_c)

# STEP 4. Structure with large pores. Pore B ________________________________________________________________________
# Load data
with open('../zif7/__tmp/atoms_zif7_lp.pickle', 'rb') as handle:
    atoms_zif7_lp = pickle.load(handle)

# Save .gro file
atoms_zif7_lp, atoms_tempo = put_tempo_in_lattice(atoms_zif7_lp, atoms_tempo, [1783, 2465])
with open('__tmp/atoms_zif7_lpB.pickle', 'wb') as handle:
    pickle.dump(atoms_zif7_lp, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('__tmp/atoms_tempo_lpB.pickle', 'wb') as handle:
    pickle.dump(atoms_tempo, handle, protocol=pickle.HIGHEST_PROTOCOL)
atoms_zif7_tempo_lpB = atoms_zif7_lp.copy()
atoms_zif7_tempo_lpB.extend(atoms_tempo.copy())
write_gro_file(atoms_zif7_tempo_lpB, 'zif7_tempo_lpB.gro', a_lp, b_lp, c_lp, alpha_lp, beta_lp, gamma_lp, bounds_a, bounds_b, bounds_c)
