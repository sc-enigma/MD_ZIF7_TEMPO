import sys
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import pickle

sys.path.append('components')
from components.IterLoad import Iterload

def compute_pos(trajectory, cnt_chunk, params):
    x = trajectory.time
    y = trajectory.xyz[:, params['ids']]
    return x, y, []

for job_name in ['zif7_tempo_lpA', 'zif7_tempo_npA', 'zif7_tempo_lpB']:
    traj = f'/media/sc_enigma/_data/Projects/MD_ZIF7/production/{job_name}/prod/traj_comp.xtc'
    top = f'/media/sc_enigma/_data/Projects/MD_ZIF7/production/{job_name}/prod/confout.gro'
    output = f'/home/sc_enigma/Projects/MD_ZIF7/analysis/msd/pos_tempo_n/{job_name}'

    topology = md.load_topology(f'/media/sc_enigma/_data/Projects/MD_ZIF7/production/{job_name}/prod/confout.gro')
    params = {}
    params['ids'] = topology.select("resname == TMP and type == N")
    params['data_type'] = 'DataTPos'
    Iterload(traj, top, output, compute_pos, params)


'''with open('/home/sc_enigma/Projects/MD_ZIF7/analysis/msd/pos_tempo_n/zif7_tempo_lpA.pickle', 'rb') as handle:
    pos_data = pickle.load(handle)
plt.plot(pos_data.x, pos_data.y[:, 0, 0], label='x')
plt.plot(pos_data.x, pos_data.y[:, 0, 1], label='y')
plt.plot(pos_data.x, pos_data.y[:, 0, 2], label='z')
plt.legend()
plt.show()'''
