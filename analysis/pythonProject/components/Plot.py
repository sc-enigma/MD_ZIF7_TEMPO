import pickle
import numpy as np
import matplotlib.pyplot as plt

# Create
# p = Plot([1,2,3], [1, 2, 3])
# p.params['linewidth'] = 10
#
# Save
# p.save('/home/sc_enigma/Projects/MD_ZIF7/analysis/p')
#
# Load
# p.load('/home/sc_enigma/Projects/MD_ZIF7/analysis/p')
# p = Plot(filename='/home/sc_enigma/Projects/MD_ZIF7/analysis/p')
#
# Show
# p.show()
# plt.show()

class Plot:
    def __init__(self, x=[], y=[], params={}, filename=''):
        if filename != '':
            self.load(filename)
        else:
            self.x = np.array(x)
            self.y = np.array(y)
            self.params = params

    def save(self, filename):
        with open(filename, 'wb') as handle:
            pickle.dump([self.x, self.y, self.params], handle, protocol=pickle.HIGHEST_PROTOCOL)

    def load(self, filename):
        with open(filename, 'rb') as handle:
            [self.x, self.y, self.params] = pickle.load(handle)

    def show(self):
        marker = None
        color = None
        linewidth = None
        label = None

        if 'marker' in self.params.keys():
            marker = self.params['marker']
        if 'color' in self.params.keys():
            color = self.params['color']
        if 'linewidth' in self.params.keys():
            linewidth = self.params['linewidth']
        if 'label' in self.params.keys():
            label = self.params['label']

        plt.plot(self.x, self.y, marker=marker, color=color, linewidth=linewidth, label=label)

# Create
# pp = PlotParams(['legend'])
# pp.params['xlim'] = [0, 10]
#
# Save
# pp.save('/home/sc_enigma/Projects/MD_ZIF7/analysis/pp')
#
# Load
# pp.load('/home/sc_enigma/Projects/MD_ZIF7/analysis/pp')
# pp = PlotParams(filename='/home/sc_enigma/Projects/MD_ZIF7/analysis/pp')
#
# Apply
# pp.apply()

class PlotParams:
    def __init__(self, params={}, filename=''):
        if filename != '':
            self.load(filename)
        else:
            self.params = params
        self.apply()

    def save(self, filename):
        with open(filename, 'wb') as handle:
            pickle.dump(self.params, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def load(self, filename):
        with open(filename, 'rb') as handle:
            self.params = pickle.load(handle)

    def apply(self):
        if 'xlim' in self.params.keys():
            plt.xlim(self.params['xlim'])
        if 'ylim' in self.params.keys():
            plt.xlim(self.params['ylim'])
        if 'legend' in self.params.keys():
            plt.legend()

