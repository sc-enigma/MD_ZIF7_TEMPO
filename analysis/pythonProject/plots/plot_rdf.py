import sys
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import pickle

sys.path.append('components')
from components.Plot import Plot

for job_name in ['zif7_tempo_lpA', 'zif7_tempo_npA', 'zif7_tempo_lpB']:
    with open('/home/sc_enigma/Projects/MD_ZIF7/analysis/rdf/zif7_tempo_npA/tmp_all_zif7_c5.pickle.pickle',
              'rb') as handle:
        rdf_data = pickle.load(handle)

with open('/home/sc_enigma/Projects/MD_ZIF7/analysis/rdf/zif7_tempo_npA/tmp_all_zif7_zn.pickle.pickle', 'rb') as handle:
    rdf_data = pickle.load(handle)
npA_zn = Plot(rdf_data.x, rdf_data.y / rdf_data.cnt_chunk)
npA_zn.params['label'] = 'TEMPO - Zn'
npA_zn.show()

with open('/home/sc_enigma/Projects/MD_ZIF7/analysis/rdf/zif7_tempo_npA/tmp_all_zif7_n.pickle.pickle', 'rb') as handle:
    rdf_data = pickle.load(handle)
npA_n = Plot(rdf_data.x, rdf_data.y / rdf_data.cnt_chunk)
npA_n.params['label'] = 'TEMPO - N'
npA_n.show()

with open('/home/sc_enigma/Projects/MD_ZIF7/analysis/rdf/zif7_tempo_npA/tmp_all_zif7_c2.pickle.pickle', 'rb') as handle:
    rdf_data = pickle.load(handle)
npA_c2 = Plot(rdf_data.x, rdf_data.y / rdf_data.cnt_chunk)
npA_c2.params['label'] = 'TEMPO - C2'
npA_c2.show()

with open('/home/sc_enigma/Projects/MD_ZIF7/analysis/rdf/zif7_tempo_npA/tmp_all_zif7_c5.pickle.pickle', 'rb') as handle:
    rdf_data = pickle.load(handle)
npA_c5 = Plot(rdf_data.x, rdf_data.y / rdf_data.cnt_chunk)
npA_c5.params['label'] = 'TEMPO - C5'
npA_c5.show()

plt.legend()
plt.xlim([0, 0.6])
plt.xlabel('r, nm')
plt.ylabel('g(r)')
plt.show()