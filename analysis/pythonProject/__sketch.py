import sys
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

sys.path.append('/home/sc_enigma/Projects/MD_ZIF7/analysis/')
data_path = '/home/sc_enigma/Projects/MD_ZIF7/analysis/'
top = md.load_topology(f'{data_path}confout.gro')

# TEMPO POSITION
'''ref_atoms = top.select("name == O1")
xyz_ref = np.empty((0, 1, 3))

for chunk in md.iterload(f'{data_path}traj_comp.xtc', top=f'{data_path}confout.gro', chunk=1000, stride=10):
    xyz_ref = np.append(xyz_ref, np.array(chunk.xyz[:, ref_atoms]), axis=0)

plt.plot(np.arange(len(xyz_ref)), xyz_ref[:,0,0])
plt.plot(np.arange(len(xyz_ref)), xyz_ref[:,0,1])
plt.plot(np.arange(len(xyz_ref)), xyz_ref[:,0,2])
plt.show()'''

# LINKER ANGLE
'''def norma(vec):
    return pow(sum(np.power(vec, 2)), 0.5)
def normalize_sgn(vec):
    length = pow(sum(np.power(vec, 2)), 0.5)
    if length > 0.0:
        vec /= length
        if vec[2] < 0.0:
            vec *= -1.0
    return vec

def angle_between(vec_a, vec_b):
    product = np.dot(vec_a, vec_b) / norma(vec_a) / norma(vec_b)
    return np.arccos(product)

angle_stats = []
zn_atoms = [[2275, 2235], [2294, 2190], [2255, 2209]]
for chunk in md.iterload(f'{data_path}traj_comp.xtc', top=f'{data_path}confout.gro', chunk=100, stride=10):
    zn_vectors = [\
        chunk.xyz[:, zn_atoms[0][0]] - chunk.xyz[:, zn_atoms[0][1]],\
        chunk.xyz[:, zn_atoms[1][0]] - chunk.xyz[:, zn_atoms[1][1]],\
        chunk.xyz[:, zn_atoms[2][0]] - chunk.xyz[:, zn_atoms[2][1]]]

    benz_vector = 0.5 * (chunk.xyz[:,2267] + chunk.xyz[:,2289]) - \
                  0.5 * (chunk.xyz[:,2274] + chunk.xyz[:,2292])

    for t_cnt in range(np.shape(zn_vectors)[1]):
        vec1 = zn_vectors[0][t_cnt]
        vec2 = zn_vectors[1][t_cnt]
        vec3 = zn_vectors[2][t_cnt]

        norm1 = normalize_sgn(np.cross(vec1, vec2))
        norm2 = normalize_sgn(np.cross(vec2, vec3))
        norm3 = normalize_sgn(np.cross(vec3, vec1))
        norm = normalize_sgn(norm1 + norm2 + norm3)
        dir = benz_vector[t_cnt]

        angle = angle_between(norm, dir)
        angle_stats.append(angle)

hist_y, hist_x = np.histogram(angle_stats, bins=50, range=[0, np.pi])
hist_y = np.append(hist_y, 0)
plt.plot(hist_x, hist_y)
plt.xlim([0, np.pi])
plt.show()'''

# WINDOW SIZE
def calc_size(vecs): # shape = (6, 3)
    return 0.0

size_stats = []
h_atoms = [2270, 2203, 2185, 2229, 2250, 2288]
for chunk in md.iterload(f'{data_path}traj_comp.xtc', top=f'{data_path}confout.gro', chunk=100, stride=10):
    h_xyz = [\
        chunk.xyz[:, h_atoms[0]],
        chunk.xyz[:, h_atoms[1]],
        chunk.xyz[:, h_atoms[2]],
        chunk.xyz[:, h_atoms[3]],
        chunk.xyz[:, h_atoms[4]],
        chunk.xyz[:, h_atoms[5]]
        ]
    for t_cnt in range(np.shape(h_xyz)[1]):
        size_stats.append(calc_size(
            [h_xyz[0][t_cnt], h_xyz[1][t_cnt], h_xyz[2][t_cnt],
            h_xyz[3][t_cnt], h_xyz[4][t_cnt], h_xyz[5][t_cnt]]))

    break

hist_y, hist_x = np.histogram(size_stats, bins=50)
hist_y = np.append(hist_y, 0)
plt.plot(hist_x, hist_y)
plt.show()