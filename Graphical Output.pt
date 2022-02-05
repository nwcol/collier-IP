#this is the script I've been using to produce plots with matplotlib. It and the main script are here as .py files because .ipynb files don't appear to 
#show up properly in github.

import statistics
import matplotlib.pyplot as plt
import numpy as np

def find_summed_variances():
    variances = []
    for gens in manygens:
        gens_variances = []
        for gen in gens:
            ancestries = []
            for ind in gen:
                ancestries.append(ind[1])
            variance = statistics.variance(ancestries)
            gens_variances.append(variance)
        variances.append(gens_variances)
    summed_variances = []
    for x in range(len(variances[0])):
        summed_variances.append(0)
    for gens in variances:
        for x in range(len(gens)):
            summed_variances[x] = summed_variances[x] + gens[x]
    for x in range(len(summed_variances)):
        summed_variances[x] = summed_variances[x] / len(manygens)
    return(summed_variances)

def find_ideal_variances():
    ideal_variances = [0.25]
    if N_const == 0:
        for x in range(1,len(manygens[0]) + 1):
            ideal_variance_x = ((S1f_fraction * (1 - S1f_fraction)) + (S1m_fraction * (1 - S1m_fraction) + (2 * c11_0))) / (2 ** (x + 1))
            ideal_variances.append(ideal_variance_x)      
    return(ideal_variances)
            
summed_variances = find_summed_variances()
ideal_variances = find_ideal_variances()

x = range(len(summed_variances))
y1 = summed_variances
y2 = ideal_variances

fig, ax = plt.subplots()
ax.plot(x, y1, y2)
ax.set_ylim(0, .25)
ax.set_xlim(0, len(x))
plt.xticks(np.arange(min(x), max(x)+1, 5.0))
leg = ax.legend(['Variance', 'Ideal Variance'])

ax.set(xlabel = 'generations', ylabel = 'variance of S1 ancestry',
      title = 'S1 Variance with N = {}, c11,0 = {}, c11 = {}, migration = {}'.format(N, c11_0, c11, N_const))
