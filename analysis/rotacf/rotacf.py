import numpy as np
import matplotlib.pyplot as plt

def read_xvg(file_name):
    file = open(file_name, 'r')
    lines = [line.replace('\n', '') for line in file]
    lines = [line for line in lines if len(line) != 0 \
             if line[0] != '@' and line[0] != '#' and line[0] != '&']
    file.close()

    data = np.transpose(np.array([[float(val) for val in line.split()] for line in lines]))
    return data[0], data[1]

def calculate_tau(x, y):
    for val_idx in range(len(x)):
        if y[val_idx] < y[0] / np.pi:
            return x[val_idx] - x[0]
    print((x[-1] - x[0]) / np.log(y[0] / y[-1]))
    return None

x, y = read_xvg('C:\\Users\\sc_enigma\\Downloads\\rotacf\\N_O\\zif7_tempo_npA.xvg')
first_idx = int(len(x) * 0.2)
x, y = x[first_idx:], y[first_idx:]
print(calculate_tau(x, y))
plt.plot(x, y)
plt.show()