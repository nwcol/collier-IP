import sys
tractspath = "C:/Genetics work/tracts-python3"
sys.path.append(tractspath)
import tracts
import pandas as pd
import tkinter as Tk
import numpy
import matplotlib.pyplot as plt
import os
import pylab
import scipy

def get_source_prop(population):
    """returns mean S1 ancestry"""
    prop_list = pop.get_means(["1"])
    S1_prop = numpy.mean(prop_list)
    return(S1_prop)
    
def init_pop(input_dir):
    file_names_all = os.listdir(input_dir)
    file_names = [file for file in file_names_all if file.split('.')[-1] == "bed"]
    names = list(set([file.split('_')[0] for file in file_names]))
    chroms = ['1']
    pop = tracts.population(names = names, fname = (input_dir, "", ".bed" ), selectchrom = chroms)
    (bins, data)=pop.get_global_tractlengths(npts=50)
    labels = ['1', '2']
    data = [data[poplab] for poplab in labels]
    return(pop, bins, data, labels)

pop, bins, data, labels = init_pop("chrom_output/")
#pop.plot_chromosome(0, {'1':'red', '2':'blue'})
#pop.plot_global_tractlengths({'1':'red', "2":"blue"})
######################################################################################################

def get_source_prop(pop):
    """returns mean S1 ancestry"""
    prop_list = pop.get_means(["1"])
    s1_prop = numpy.mean(prop_list)
    return(s1_prop)

def test_func(params, fracs):
    """a model to figure out how models work"""
    t_start = params[0] * 100
    if t_start < 0:
        return("t_start went below 0!")
    mig  = numpy.zeros((int(numpy.ceil(t_start)) + 1, 2))
    mig[-1,:] = fracs
    gen = int(numpy.ceil(t_start)) + 1
    scale = gen - t_start - 1
    mig[-2,:] = scale * fracs
    return(mig)

def test_func_constraint(params, fracs):
    """constraint function"""
    ret = 1
    (t_start) = params
    t_start = t_start * 100
    (s1_init, s2_init) = fracs
    ret = min(1, 1 - s1_init) #migration proportion must be between 0 and 1
    ret = min(ret, s1_init)
    # generate the migration matrix and test for possible issues
    func = test_func #specify the model
    mig = func(params, fracs) #get the migration matrix
    totmig = mig.sum(axis=1) #calculate the migration rate per generation
    ret = min(ret, -abs(totmig[-1] - 1) + 1e-8) #first generation migration must sum up to 1
    ret = min(ret, -totmig[0], -totmig[1]) #no migrations are allowed in the first two generations
    #migration at any given generation cannot be greater than 1
    ret = min(ret, 10*min(1-totmig), 10*min(totmig))
    #start time must be at least two generations ago
    ret = min(ret, t_start - 0.02)
    # print some diagnistics (facultative)
    if abs(totmig[-1]-1) > 1e-8:
        print(mig)
        print("founding migration should sum up to 1. Now:")
    if totmig[0] > 1e-10:
        print("migrants at last generation should be removed from sample!")
    if totmig[1] > 1e-10:
        print("migrants at penultimate generation should be removed from sample!")
    if ((totmig > 1).any() or (mig < 0).any()):
        print("migration rates should be between 0 and 1")
    return ret
    
#####
t_start = 0.02
params = numpy.array([t_start])
s1_init = get_source_prop(pop)
fracs = numpy.array([s1_init, 1 - s1_init])
#####

demog = tracts.demographic_model(test_func(params, fracs))
nsamp = pop.nind
Ls = pop.Ls
xopt = tracts.optimize_cob_fracs2(
    params, 
    bins, 
    Ls, 
    data, 
    nsamp, 
    test_func, 
    fracs, 
    outofbounds_fun = test_func_constraint,
    cutoff = 1, 
    epsilon = 0.01)  

out = test_func(xopt, fracs)
print(out)
print("generations:", len(out) - 1)
