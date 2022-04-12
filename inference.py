import numpy as np
import os
import sys
tractspath = "C:/Genetics work/tracts-python3"
sys.path.append(tractspath)
import tracts


class TwoParamInference:
    
    def __init__(self, input_dir):
        self.input_dir = input_dir
        
    def make_inference(self):
        pop, bins, data, labels = self.init_pop(self.input_dir)
        mig, xopt = self.optimize(pop, bins, data, labels)
        
        (t_start, s1_init) = xopt
        return(mig, round(t_start * 100, 3), round(s1_init, 3), round(1 - s1_init, 3))

    def init_pop(self, input_dir):
        file_names_all = os.listdir(input_dir)
        file_names = [file for file in file_names_all if file.split('.')[-1] == "bed"]
        names = list(set([file.split('_')[0] for file in file_names]))
        chroms = ['1']
        pop = tracts.population(names = names, fname = (input_dir, "", ".bed" ), selectchrom = chroms)
        (bins, data)=pop.get_global_tractlengths(npts=50)
        labels = ['1', '2']
        data = [data[poplab] for poplab in labels]
        return(pop, bins, data, labels)
    
    def get_source_prop(self, pop):
        """returns mean S1 ancestry"""
        prop_list = pop.get_means(["1"])
        s1_prop = np.mean(prop_list)
        return(s1_prop)

    def two_param_model(self, params):
        """a model which optimizes two parameters:
            generations since initial admixture
            initial s1 (and by extension initial s2, since s2 = 1 - s1)
        """
        t_start = params[0] * 100
        gen = int(np.ceil(t_start)) + 1
        init_props = np.array([params[1], 1 - params[1]])
        mig  = np.zeros([gen, 2])
        mig[-1,:] = init_props
        scale = gen - t_start - 1
        mig[-2,:] = scale * init_props
        return(mig)

    def two_param_constraint(self, params):
        """constraint function"""
        ret = 1
        (t_start, s1_init) = params
        mig = self.two_param_model(params)
        totmig = mig.sum(axis=1)
        ret = min(ret, -abs(totmig[-1]-1)+1e-8)
        ret = min(ret, -totmig[0], -totmig[1])
        ret = min(ret, 10*min(1-totmig), 10*min(totmig))
        ret = min(ret, t_start-.02)
        if abs(totmig[-1]-1) > 1e-8:
            print(mig)
            print("founding migration should sum up to 1. Now:")
        if totmig[0] > 1e-10:
            print("migrants at last generation should be removed from sample!")
        if totmig[1] > 1e-10:
            print("migrants at penultimate generation should be removed from sample!")
        if ((totmig > 1).any() or (mig < 0).any()):
            print("migration rates should be between 0 and 1")
        return(ret)

    def optimize(self, pop, bins, data, labels):
        t_start = 0.03
        s1_init = self.get_source_prop(pop)
        params = np.array([t_start, s1_init])
        demog = tracts.demographic_model(self.two_param_model(params))
        nsamp = pop.nind
        Ls = pop.Ls
        xopt = tracts.optimize_cob(
            p0 = params, 
            bins = bins, 
            Ls = Ls, 
            data = data, 
            nsamp = nsamp, 
            model_func = self.two_param_model,  
            outofbounds_fun = self.two_param_constraint,
            verbose = 0)  
        mig = self.two_param_model(xopt)
        return(mig, xopt)


class FourParamInference:
    
    def __init__(self, input_dir):
        self.input_dir = input_dir
        
    def make_inference(self):
        pop, bins, data, labels = self.init_pop(self.input_dir)
        mig, xopt = self.optimize(pop, bins, data, labels)
        
        (t_start, s1_init, s1_const, s2_const) = xopt
        gen = t_start * 100
        return(mig, "gen:", gen, "s1_init:", s1_init, "s2_init:", 1 - s1_init,
                "s1_const:", s1_const, "s2_const:", s2_const)

    def init_pop(self, input_dir):
        file_names_all = os.listdir(input_dir)
        file_names = [file for file in file_names_all if file.split('.')[-1] == "bed"]
        names = list(set([file.split('_')[0] for file in file_names]))
        chroms = ['1']
        pop = tracts.population(names = names, fname = (input_dir, "", ".bed" ), selectchrom = chroms)
        (bins, data)=pop.get_global_tractlengths(npts=50)
        labels = ['1', '2']
        data = [data[poplab] for poplab in labels]
        return(pop, bins, data, labels)
    
    def get_source_prop(self, pop):
        """returns mean S1 ancestry"""
        prop_list = pop.get_means(["1"])
        s1_prop = np.mean(prop_list)
        return(s1_prop)

    def four_param_model(self, params):
        """a model which optimizes five parameters:
            generations since initial admixture
            initial S1/S2 contribution
            constant (migrant) S1/S2 contribution
        """
        t_start = params[0] * 100
        init_props = np.array([params[1], 1 - params[1]])
        if t_start < 0:
            return("t_start went below 0!")
        mig  = np.zeros((int(np.ceil(t_start)) + 1, 2))
        mig[:,:1] = params[2]
        mig[:,-1] = params[3]
        mig[-1,:] = init_props
        gen = int(np.ceil(t_start)) + 1
        scale = gen - t_start - 1
        mig[-2,:] = scale * init_props
        mig[:2,:] = np.array([0, 0])
        return(mig)

    def four_param_constraint(self, params):
        """constraint function"""
        ret = 1
        (t_start, s1_init, s1_const, s2_const) = params
        
        mig = self.four_param_model(params) #get the migration matrix
        totmig = mig.sum(axis=1) #calculate the migration rate per generation

        ret = min(ret, -abs(totmig[-1]-1)+1e-8) #first generation migration must sum up to 1
        ret = min(ret, -totmig[0], -totmig[1]) #no migrations are allowed in the first two generations

        #migration at any given generation cannot be greater than 1
        ret = min(ret, 10*min(1-totmig), 10*min(totmig))

        #start time must be at least two generations ago
        ret = min(ret, t_start-.02)

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

        return(ret)

    def optimize(self, pop, bins, data, labels):
        t_start = 0.02
        s1_init = 0.50
        s1_const = 0.1
        s2_const = 0.1
        params = np.array([t_start, s1_init, s1_const, s2_const])

        demog = tracts.demographic_model(self.four_param_model(params))
        nsamp = pop.nind
        Ls = pop.Ls
        xopt = tracts.optimize_cob(
            p0 = params, 
            bins = bins, 
            Ls = Ls, 
            data = data, 
            nsamp = nsamp, 
            model_func = self.four_param_model,  
            outofbounds_fun = self.four_param_constraint)  

        mig = self.four_param_model(xopt)
        return(mig, xopt)