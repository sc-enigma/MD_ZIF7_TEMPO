import sys
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import pickle

sys.path.append('components')
from components.IterLoad import Iterload

def compute_rdf(trajectory, cnt_chunk, params):
    x, y = md.compute_rdf(trajectory, params['pairs'], r_range=[0.0, 2.0], bin_width=0.005, n_bins=None, periodic=True, opt=True)
    return x, y, []

for job_name in ['zif7_tempo_lpB']: # 'zif7_tempo_lpA', 'zif7_tempo_npA'
    traj = f'/media/sc_enigma/_data/Projects/MD_ZIF7/production/{job_name}/prod/traj_comp.xtc'
    top = f'/media/sc_enigma/_data/Projects/MD_ZIF7/production/{job_name}/prod/confout.gro'

    # dict val = [x, y, cnt]
    pairs = {}
    topology = md.load_topology(f'/media/sc_enigma/_data/Projects/MD_ZIF7/production/{job_name}/prod/confout.gro')
    pairs['tmp_all_zif7_zn'] = topology.select_pairs("resname == TMP", "name == Zn")
    pairs['tmp_all_zif7_n'] = topology.select_pairs("resname == TMP", "name == N")
    pairs['tmp_all_zif7_c2'] = topology.select_pairs("resname == TMP", "name == C2")
    pairs['tmp_all_zif7_c5'] = topology.select_pairs("resname == TMP", "name == C5")
    pairs['tmp_all_zif7_c6'] = topology.select_pairs("resname == TMP", "name == C6")
    pairs['tmp_n_zif7_zn'] = topology.select_pairs("resname == TMP and type == N", "name == Zn")
    pairs['tmp_n_zif7_n'] = topology.select_pairs("resname == TMP and type == N", "name == N")
    pairs['tmp_n_zif7_c2'] = topology.select_pairs("resname == TMP and type == N", "name == C2")
    pairs['tmp_n_zif7_c5'] = topology.select_pairs("resname == TMP and type == N", "name == C5")
    pairs['tmp_n_zif7_c6'] = topology.select_pairs("resname == TMP and type == N", "name == C6")

    for key in pairs.keys():
        print(job_name, key)
        output = f'/home/sc_enigma/Projects/MD_ZIF7/analysis/rdf/{job_name}/{key}'
        params = {}
        params['pairs'] = pairs[key]
        Iterload(traj, top, output, compute_rdf, params)

'''
with open('/home/sc_enigma/Projects/MD_ZIF7/analysis/rdf/zif7_tempo_lpA/zif7_tempo_lpA.pickle.pickle', 'rb') as handle:
    rdf_data = pickle.load(handle)
plt.plot(rdf_data.x, rdf_data.y / rdf_data.cnt_chunk)
plt.show()
'''