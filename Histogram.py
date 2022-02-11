#Creates a series of histograms of gen0 to gen6 which are very helpful in visualising and interpreting the patterns of admixture

import numpy as np

import matplotlib.pyplot as plt

fig, ax = plt.subplots(7, 1, figsize = (5.9, 46))

gens = manygens[0]
for x in range(0,7):
    ancestries = []
    for y in gens[x]:
        ancestries.append(y[1]) 
        
    variance = round(np.var(ancestries), 5)
    ax[x].hist(ancestries, bins = 100, range = (0, 1), density = True, color = 'red')
    ax[x].set(title = ('gen = {}, variance = {}'.format(x, variance)))
    ax[x].set_xlim(-0.05, 1.05)
    ax[x].set_ylim(0, 100)

print('done')
