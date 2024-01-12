import sys
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import pickle
from numba import njit

sys.path.append('components')
from components.Plot import Plot

@njit
def average_pos(y, delta=500):
    y_avg = np.zeros(len(y))
    for i in range(len(y)):
        j_min = max(0, i - delta)
        j_max = min(len(y) - 1, i + delta)
        cnt = 0
        for j in range(j_min, j_max):
            y_avg[i] += y[j]
            cnt += 1
        y_avg[i] /= cnt
    return y_avg

log_file = open('/home/sc_enigma/Projects/MD_ZIF7/analysis/msd/pos_tempo_n/log.txt', 'w')
for job_name in ['zif7_tempo_lpA', 'zif7_tempo_npA', 'zif7_tempo_lpB']:
    with open(f'/home/sc_enigma/Projects/MD_ZIF7/analysis/msd/pos_tempo_n/{job_name}.pickle', 'rb') as handle:
        pos_data = pickle.load(handle)

    t = pos_data.x
    x = pos_data.y[:, 0, 0]
    y = pos_data.y[:, 0, 1]
    z = pos_data.y[:, 0, 2]
    x_avg = average_pos(x)
    y_avg = average_pos(y)
    z_avg = average_pos(z)

    '''plt.plot(t, x, label='x')
    plt.plot(t, y, label='y')
    plt.plot(t, z, label='z')
    plt.plot(t, x_avg, label='x_avg')
    plt.plot(t, y_avg, label='y_avg')
    plt.plot(t, z_avg, label='z_avg')
    plt.legend()
    plt.show()'''

    msd = np.sqrt(sum(np.power(x_avg - x, 2) + np.power(y_avg - y, 2) + np.power(z_avg - z, 2)) / len(x))
    print(f'MSD {job_name} = {round(msd, 3)} nm')
    log_file.write(f'MSD {job_name} = {round(msd, 3)} nm\n')

log_file.close()