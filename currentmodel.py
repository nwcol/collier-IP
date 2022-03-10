import numpy as np
import tskit
import msprime
import matplotlib.pyplot as plt
import random
import pandas as pd
import copy
from IPython.display import SVG

class Model:  
    
    def __init__(
        self, 
        g,
        description = None, 
        runs = 1,  
        S1f_init = 0, 
        S1m_init = 0, 
        S2f_init = 0, 
        S2m_init = 0, 
        N_init = None,
        S1f_const = 0, 
        S1m_const = 0, 
        S2f_const = 0, 
        S2m_const = 0,
        N_const =None,
        c11_0 = 0, 
        c22_0 = 0, 
        c11 = 0, 
        c1h = 0, 
        c12 = 0,
        ch1 = 0, 
        chh = 0, 
        ch2 = 0,
        c21 = 0, 
        c2h = 0, 
        c22 = 0):
        """Initialize the model.
            arguments
                description (string) a short description of the 
                    parameters of interest (optional)
                runs (int) the number of trials to perform (optional)
                g (int) the number of generations to simulate
                Six_init (int) the initial population by population 
                    of origin and sex
                Six_const (int) the number of migrants per generation 
                    by population of origin and sex (optional)
                cij_0 (float) the assortative constant for the first 
                    generation
                cij (float) the assortative constnat for subsequent 
                    generations
        """
        self.description = description
    
        self.runs = runs
        self.g = g
        
        if N_init == None:
            self.S1f_init = S1f_init
            self.S1m_init = S1m_init
            self.S2f_init = S2m_init
            self.S2m_init = S2f_init
        else:
            self.S1f_init = N_init // 4
            self.S1m_init = N_init // 4
            self.S2f_init = N_init // 4
            self.S2m_init = N_init // 4
            
        if N_const == None:
            self.S1f_const = S1f_const
            self.S1m_const = S1m_const
            self.S2f_const = S2f_const
            self.S2m_const = S2m_const
        else: 
            self.S1f_const = N_const // 4
            self.S1m_const = N_const // 4
            self.S2f_const = N_const // 4
            self.S2m_const = N_const // 4
        
        self.c11_0 = c11_0
        self.c12_0 = -c11_0
        self.c22_0 = c22_0
        self.c21_0 = -c22_0
        
        self.c11 = c11
        self.c1h = c1h
        self.c12 = c12
        self.ch1 = ch1
        self.chh = chh
        self.ch2 = ch2
        self.c21 = c21
        self.c2h = c2h
        self.c22 = c22

        if c11 + c1h + c12 != 0 or\
        ch1 + chh + ch2 != 0 or\
        c21 + c2h + c22 != 0 or\
        c11 + ch1 + c21 != 0 or\
        c1h + chh + c2h != 0 or\
        c21 + c2h + c22 != 0:
            raise Exception('c values do not sum to zero!')
        
        self.N_init = (self.S1f_init + self.S1m_init
            + self.S2f_init + self.S2m_init)
        self.N_const = (self.S1f_const + self.S1m_const 
            + self.S2f_const + self.S2m_const)
        self.N_f_init = self.S1f_init + self.S2f_init
        self.N_m_init = self.S1m_init + self.S2m_init 
        
    def set_c(self):
        """calculates all c values using the input cs"""
        pass
        
    def new_population(self):
        """ generates an initial generation 0 based on initial 
            parameters
        """
        gen_0 = []
        for x in range(self.S1f_init):
            gen_0.append([0, 1, 1])
        for x in range(self.S1m_init):
            gen_0.append([1, 1, 1])
        for x in range(self.S2f_init):
            gen_0.append([0, 0, 2])
        for x in range(self.S2m_init):
            gen_0.append([1, 0, 2])
        for x in range(len(gen_0)):
            gen_0[x].append(x)   
        self.id_shift = len(gen_0)
            
        return(gen_0)
    
    def const_mig(self, gen_g_minus1):
        """appends constant migrations to gen_g_minus1"""
        for x in range(self.S1f_const):
            gen_g_minus1.append([0, 1, 1, -1])
        for x in range(self.S1m_const):
            gen_g_minus1.append([1, 1, 1, -1])
        for x in range(self.S2f_const):
            gen_g_minus1.append([0, 0, 2, -2])
        for x in range(self.S2m_const):
            gen_g_minus1.append([1, 0, 2, -2])
        return(gen_g_minus1)
    
    def find_probs_gen1(self, gen_0):
        """returns probabilities of matings for gen 1"""
        S1f, S2f, S1m, S2m = 0, 0, 0, 0
        Nf, Nm = 0, 0
        for ind in gen_0:
            if ind[0] == 0:
                Nf = Nf + 1
                if ind[2] == 1:
                    S1f = S1f + 1
                elif ind[2] == 2:
                    S2f = S2f + 1
            elif ind[0] == 1:
                Nm = Nm + 1
                if ind[2] == 1:
                    S1m = S1m + 1
                elif ind[2] == 2:
                    S2m = S2m + 1

        s1f = S1f/Nf
        s2f = S2f/Nf
        s1m = S1m/Nm 
        s2m = S2m/Nm

        P11 = s1f * s1m + self.c11_0
        P12 = s1f * s2m + self.c12_0
        P21 = s2f * s1m + self.c21_0
        P22 = s2f * s2m + self.c22_0
        P1h, Ph1, Phh, Ph2, P2h = 0, 0, 0, 0, 0
        prob_list = [P11, P1h, P12, Ph1, Phh, Ph2, P21, P2h, P22]
        
        return(prob_list)

    def find_probs(self, gen_g_minus1):
        """ returns probabilities of matings based on composition 
            of gen g-1
        """
        S1f, Hf, S2f, S1m, Hm, S2m = 0, 0, 0, 0, 0, 0
        Nf, Nm = 0, 0
        for ind in gen_g_minus1:
            if ind[0] == 0:
                Nf = Nf + 1
                if ind[2] == 1:
                    S1f = S1f + 1
                elif ind[2] == 0:
                    Hf = Hf + 1
                elif ind[2] == 2:
                    S2f = S2f + 1
            elif ind[0] == 1:
                Nm = Nm + 1
                if ind[2] == 1:
                    S1m = S1m + 1
                elif ind[2] == 0:
                    Hm = Hm + 1
                elif ind[2] == 2:
                    S2m = S2m + 1

        s1f = S1f/Nf
        hf = Hf/Nf
        s2f = S2f/Nf
        s1m = S1m/Nm
        hm = Hm/Nm
        s2m = S2m/Nm

        P11 = s1f * s1m + self.c11
        P1h = s1f * hm + self.c1h
        P12 = s1f * s2m + self.c12
        Ph1 = hf * s1m + self.ch1
        Phh = hf * hm + self.chh
        Ph2 = hf * s2m + self.ch2
        P21 = s2f * s1m + self.c21
        P2h = s2f * hm + self.c2h
        P22 = s2f * s2m + self.c22
        prob_list = [P11, P1h, P12, Ph1, Phh, Ph2, P21, P2h, P22]
     
        return(prob_list)

    def mate(self, gen_g_minus1, is_gen_1, g):
        """creates a new generation"""
        if is_gen_1 == True:
            prob_list = self.find_probs_gen1(gen_g_minus1)
        else:
            prob_list = self.find_probs(gen_g_minus1)
        S1f_inds, S1m_inds, Hf_inds = [], [], []
        Hm_inds, S2f_inds, S2m_inds = [], [], []
        for ind in gen_g_minus1:
            if ind[0] == 0:
                if ind[2] == 1:
                    S1f_inds.append(ind)
                elif ind[2] == 0:
                    Hf_inds.append(ind)
                elif ind[2] == 2:
                    S2f_inds.append(ind)     
            elif ind[0] == 1:
                if ind[2] == 1:
                    S1m_inds.append(ind)
                elif ind[2] == 0:
                    Hm_inds.append(ind)
                elif ind[2] == 2:
                    S2m_inds.append(ind)

        
        matings = ['11', '1h', '12', 'h1', 'hh', 'h2', '21', '2h', '22']
        gen_g = []
        gen_g_parents = []
        for x in range(self.N_init):
            
            [mating] = random.choices(matings, weights = prob_list)
            
            if mating == '11':
                female_parent = random.choice(S1f_inds)
                male_parent = random.choice(S1m_inds)
                #parent index in ancestries and offspring index should match.
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] 
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '1h':
                female_parent = random.choice(S1f_inds)
                male_parent = random.choice(Hm_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] 
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '12':
                female_parent = random.choice(S1f_inds)
                male_parent = random.choice(S2m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] 
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == 'h1':
                female_parent = random.choice(Hf_inds)
                male_parent = random.choice(S1m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] 
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == 'hh':
                female_parent = random.choice(Hf_inds)
                male_parent = random.choice(Hm_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] 
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring)


            elif mating == 'h2':
                female_parent = random.choice(Hf_inds)
                male_parent = random.choice(S2m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] 
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '21':
                female_parent = random.choice(S2f_inds)
                male_parent = random.choice(S1m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] 
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '2h':
                female_parent = random.choice(S2f_inds)
                male_parent = random.choice(Hm_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1]
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '22':
                female_parent = random.choice(S2f_inds)
                male_parent = random.choice(S2m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] 
                    + male_parent[1]) / 2, 0]
                gen_g.append(offspring) 
        
        start_id = self.N_init * g
        for x in range(len(gen_g)):
            gen_g[x].append(start_id + x)

        return(gen_g, gen_g_parents)

    def tabulate_pedigree_gen_0(self, gen_0):
        """tabulates parentage for gen 0 individuals. all parents in 
            gen 0 are -1 or -2
        """
        ids = []
        ancestries = []
        pops = []
        female_parents = []
        male_parents = []

        for index in range(len(gen_0)):
            ind = gen_0[index]
            pops.append(ind[2])
            ids.append(ind[3])
            if ind[2] == 1:
                ancestries.append(1)
                female_parents.append(-1)
                male_parents.append(-1)
            elif ind[2] == 2:
                ancestries.append(0)
                female_parents.append(-2)
                male_parents.append(-2)

        gen_0_pedigree_table = {'g' : [0] * len(gen_0),
                              'id' : range(len(gen_0)),
                              'ancestry' : ancestries,
                              'pop' : pops,
                              'female_parent' : female_parents,
                              'male_parent' : male_parents   }
        gen_0_pedigree_df = pd.DataFrame(gen_0_pedigree_table)

        return(gen_0_pedigree_df)

    def tabulate_pedigree(self, gen_g, gen_g_parents, x):
        """ tabulates parentage for gen g individuals. x is 
            generation number
        """
        ids = []
        ancestries = []
        pops = []
        female_parents = []
        male_parents = []

        for index in range(len(gen_g)):
            parents = gen_g_parents[index]
            female_parent = parents[0]
            male_parent = parents[1]
            female_parents.append(female_parent[3]) 
            male_parents.append(male_parent[3])

            ind = gen_g[index]
            pops.append(ind[2])
            ancestries.append(ind[1])
            ids.append(ind[3])

        gen_g_pedigree_table = {'g' : [x] * len(gen_g),
                              'id' : ids,
                              'ancestry' : ancestries,
                              'pop' : pops,
                              'female_parent' : female_parents,
                              'male_parent' : male_parents   }
        gen_g_pedigree_df = pd.DataFrame(gen_g_pedigree_table)

        return(gen_g_pedigree_df)

    def make_gens(self):
        """produces a list of generations"""        
        gens = []
        internal_pedigree_list = []

        gen_0 = self.new_population()
        gens.append(gen_0) 
        gen_g_minus1 = copy.deepcopy(gen_0)
        gen_0_pedigree_df = self.tabulate_pedigree_gen_0(gen_0)
        internal_pedigree_list.append(gen_0_pedigree_df)

        for x in range(1, self.g + 1):
            if x == 1:
                gen_g, gen_g_parents = self.mate(gen_g_minus1, True, x)
                gens.append(gen_g)
            else:
                del(gen_g_minus1[(self.N_init - self.N_const):])
                gen_g_minus1 = self.const_mig(gen_g_minus1)
                gen_g, gen_g_parents = self.mate(gen_g_minus1, False, x)
                gens.append(gen_g)

            gen_g_pedigree_df = self.tabulate_pedigree(gen_g, 
                gen_g_parents, x)
            internal_pedigree_list.append(gen_g_pedigree_df)
            gen_g_minus1 = copy.deepcopy(gen_g)

        pedigree_df = pd.concat(internal_pedigree_list)
        pedigree_df = pedigree_df.reset_index(drop = True) 

        return(gens, pedigree_df)

    def execute_model(self):
        """Integrates functions to perform a series of trials"""
        gens_list = []
        pedigree_df_list = []

        for x in range(self.runs):
            gens, pedigree_df = self.make_gens()
            gens_list.append(gens)
            pedigree_df_list.append(pedigree_df)
        
        self.gens_list = gens_list
        self.pedigree_df_list = pedigree_df_list 
        
    def var_plot(self):
        
        var = []
        testgens = self.gens_list[0]
        for x in range(len(testgens)):
            gen_x_ancestries = []
            for gens in self.gens_list:
                current_gen = gens[x]
                for ind in current_gen:
                    gen_x_ancestries.append(ind[1])
            gen_x_var = np.var(gen_x_ancestries)
            var.append(gen_x_var)

        var[0] = None

        x = range(len(var))

        varx.plot(x, var)
        varx.set_ylim(0, .25)
        varx.set_xlim(0, len(x))
        plt.xticks(np.arange(min(x), max(x) + 1, 5))
    
    def histogram(self):
        
        fig, ax = plt.subplots(self.g, 1, figsize = (5.9, 7 * self.g))

        gens = self.gens_list[0]
        for x in range(0, self.g):
            ancestries = []
            for y in gens[x]:
                ancestries.append(y[1]) 

            variance = round(np.var(ancestries), 5)
            
            ax[x].hist(
                ancestries, 
                bins = 100, 
                range = (0, 1), 
                density = True, 
                color = 'red')
            ax[x].set(title = ('gen = {}, variance = {}'.format(x, variance)))
            ax[x].set_xlim(-0.05, 1.05)
            ax[x].set_ylim(0, 100)  
        
    def compare_predicted_var(self):
        
        var = []
        testgens = self.gens_list[0]
        for x in range(len(testgens)):  
            gen_x_ancestries = []
            for gens in self.gens_list:
                current_gen = gens[x]
                for ind in current_gen:
                    gen_x_ancestries.append(ind[1])    
            gen_x_var = np.var(gen_x_ancestries)
            var.append(gen_x_var)   
        var[0] = None
        
        s1f_0 = self.S1f_init / self.N_f_init  
        s1m_0 = self.S1m_init / self.N_m_init
        s2f_0 = self.S2f_init / self.N_f_init
        s2m_0 = self.S2m_init / self.N_m_init

        c11_0 = self.c11_0
        c12_0 = self.c12_0
        c21_0 = self.c21_0

        s1f = self.S1f_const / self.N_f_init 
        s1m = self.S1m_const / self.N_m_init 
        hf = (self.N_f_init - self.S1f_const 
            - self.S2f_const) / self.N_f_init
        hm = (self.N_m_init - self.S1m_const 
            - self.S2m_const) / self.N_m_init
        s2f = self.S2f_const / self.N_f_init
        s2m = self.S2m_const / self.N_m_init

        s1_0 = (s1f_0 + s1m_0) / 2
        s1 = (s1f + s1m) / 2
        s2 = (s2f + s2m) / 2
        h = 1 - s1 - s2 

        c11 = self.c11
        c1h = self.c1h
        c12 = self.c12
        ch1 = self.ch1
        chh = self.chh
        ch2 = self.ch2
        c21 = self.c21
        c2h = self.c2h
        c22 = self.c22  #doesn't appear in functions?

        predicted_var = []
        predicted_var.append(None)
        E_H_list = []
        E_H_list.append(None)
        
        for g in range(1, 21):
            if g == 1:
                #from eq 16
                E_H_list.append((s1f_0 * s1m_0) + c11_0 + (((s1f_0 * s2m_0) + 
                    (s2f_0 * s1m_0) + c12_0 + c21_0) / 2))
                #from eq 20
                predicted_var.append(((s1f_0 * (1 - s1f_0)) + 
                    (s1m_0 * (1 - s1m_0)) + (2 * c11_0)) / 4)
            elif g > 1:
                #from eq 17
                E_H_list.append(
                    ((s1f * s1m) + c11) +
                    (((s1f * s2m) + (s2f * s1m) + c12 + c21 + c1h + ch1) / 2) +
                    ((((s1f * hm) + (hf * s1m) + c1h + ch1) / 2) * (1 + E_H_list[g-1])) +
                    ((((hf * hm) + chh) / 2) * (2 * E_H_list[g-1])) +
                    ((((s2f * hm) + (hf * s2m) + ch2 + c2h) / 2) * E_H_list[g-1]))                     
                #from eq 21
                predicted_var.append(
                    (((s1f * (1 - s1f)) + (s1m * (1 - s1m)) + (2 * c11)) / 4) +
                    (((c1h + ch1 - (s1f * hf) - (s1m * hm)) / 2) * E_H_list[g-1]) +
                    ((((hf * (1 - hf)) + (hm * (1 - hm)) + (2 * chh)) / 4) * (E_H_list[g-1] ** 2)) +
                    (((hf + hm) / 4) * predicted_var[g-1]))
        
        x = range(len(var))
        predicted_varfig = plt.figure()
        var_predx = predicted_varfig.add_subplot(111)
        var_predx.plot(x, var, color = ('red'))
        var_predx.plot(predicted_var, color = ('black'))
        var_predx.set_ylim(0, .25)
        var_predx.set_xlim(0, len(x))
        plt.xticks(np.arange(min(x), max(x) + 1, 5))
        
        var_predx.set(
            xlabel = 'generations', 
            ylabel = 'variance of S1 ancestry', 
            title = 'Variance in Ancestry')
        var_predx.legend([self.description, 'predicted variance']) 
    
    def write_pedigree_df(self, output_filename, df_index = 0):
        
        selected_df = self.pedigree_df_list[df_index]
        selected_df.to_csv(output_filename) 
        return("file output successful")
    
    def write_pedigree_df_list(self, output_filename = None):
        """h"""
        if output_filename == None:
            output_filename = ("model df_list N" + str(self.N_init) + " m" + 
                str(self.N_const) + ".csv")
            output_filename = "modeldata/" + output_filename     
        output_file = open(output_filename, 'w')
        for pedigree_df in self.pedigree_df_list:
            pedigree_df.to_csv(output_file)
        output_file.close
        return("file output successful")
    
    def re_index(self, df):
        '''re-indexes sample_pedigree so that it can be used in msprime'''
        for row in range(len(df)):
            df.iloc[row, 0] = str(df.iloc[row, 0])
        for row in range(len(df)):
            replaced_id = df.iloc[row, 1]
            df = df.replace(to_replace = replaced_id, value = str(row))
            
        return(df)
    
    def format_df(self, df):
        '''corrects variable types for the sample pedigree df'''
        df['g'] = pd.to_numeric(df['g'])
        df['id'] = pd.to_numeric(df['id'])
        df['pop'] = pd.to_numeric(df['pop'])
        df['female_parent'] = pd.to_numeric(df['female_parent'])
        df['male_parent'] = pd.to_numeric(df['male_parent'])
        
        return(df)

    def sample_pedigree_df(self, sample_size, df_index = 0):
        """samples a specificed number of individuals from the last 
            generation and trims all individuals which aren't 
            ancestral to them out"""
        #create the sample list
        pedigree_df = self.pedigree_df_list[df_index]
        last_gen_list =[]
        for row in range(len(pedigree_df)):
            if pedigree_df.iloc[row, 0] == self.g:
                last_gen_list.append(pedigree_df.iloc[[row]])
        sample = random.sample(last_gen_list, sample_size)
        sample_list = sample
        
        for g in range(self.g - 1, -1, -1):
            sample_parent_ids = []
            for ind in sample:
                sample_parent_ids.append(ind.iloc[0, 4])
                sample_parent_ids.append(ind.iloc[0, 5])
            #remove duplicate parent ids
            sample_parent_ids = list(set(sample_parent_ids))
            sample_parents = []
            for parent_id in sample_parent_ids:
                if parent_id > 0:
                    sample_parents.append(pedigree_df.iloc[[parent_id]])
                else:
                    pass
            sample_list = sample_list + sample_parents
            sample = sample_parents
        
        sample_df = pd.concat(sample_list)
        sample_df = sample_df.sort_values(by = ["g", "id"])
        sample_df = self.re_index(sample_df)
        sample_df = self.format_df(sample_df)
        sample_df = sample_df.reset_index(drop = True)
        
        return(sample_df)
      
    def write_sample_df(self, sample_size, output_filename, df_index = 0):
        """write sample_df to file as .csv"""
        sample_df = self.sample_pedigree_df(sample_size, df_index)
        sample_df.to_csv(output_filename) 
        
        return("file output successful")
        
    def write_sample_df_as_txt(self, output_filename):
        """currently unused"""
        sample_df = self.sample_df
        for row in range(len(self.sample_df)):
            sample_df.iloc[row, 0] = 5 - sample_df.iloc[row, 0]
        sample_df.drop('ancestry', inplace = True, axis = 1)
        sample_df = sample_df.rename(columns={
            "g" : "time", 
            "id" :"id", 
            "pop" :"location", 
            "female_parent" : "parent0", 
            "male_parent" : "parent1"})
        sample_df = sample_df[[
            "id", 
            "location", 
            "parent0", 
            "parent1", 
            "time"]]
        sample_string = sample_df.to_string(index = False)
        output_filename = "modeldata/" + output_filename + ".txt"
        output_file = open(output_filename, "w")
        output_file.write(sample_string)
        output_file.close()
        
        return("file output successful")

    
class Sim:
   
    def __init__(self, input_filename, seq_length, recomb_rate):
        """
            arguments: 
                input_filename (string) specifies the path and file 
                    name from which the pedigree file is loaded
                seq_length (int) specifies the length of simulated 
                    chromosomes in base pairs
                recomb_rate (float) specifies the chance of 
                    recombination between neighboring base pairs       
        """
        self.input_filename = input_filename
        self.seq_length = seq_length
        self.recomb_rate = recomb_rate
        
    def make_pedigree(self):
        """ scans the pedigree table created by model.py and creates 
            a pedigree tree structure.
        """
        pedigree_df = pd.read_csv(self.input_filename)
        max_gen = pedigree_df.iloc[len(pedigree_df) - 1, 1]
        self.max_gen = max_gen
        demography = msprime.Demography()
        demography.add_population(name = "H", initial_size = 0)
        demography.add_population(name = "S1", initial_size = 0)
        demography.add_population(name = "S2", initial_size = 0)
        self.demography = demography
        this_pedigree = msprime.PedigreeBuilder(demography = self.demography)
        
        #individual 0 or nodes 0, 1 is the virtual the progenitor of all S1
        #indivudal 1 or nodes 2, 3 is the virtual progenitor of all S2
        this_pedigree.add_individual(
            time = max_gen + 1, 
            parents = [-1, -1], 
            population = "S1")
        this_pedigree.add_individual(
            time = max_gen + 1, 
            parents = [-1, -1], 
            population = "S2")
        
        for x in range(len(pedigree_df)):
            female_parent_id = pedigree_df.iloc[x, 5] + 2
            male_parent_id = pedigree_df.iloc[x, 6] + 2
            if female_parent_id == -1:
                female_parent = 0
            elif female_parent_id == -2:
                female_parent = 1
            else:
                female_parent = female_parent_id
            if male_parent_id == -1:
                male_parent = 0
            elif male_parent_id == -2:
                male_parent = 1
            else:
                male_parent = male_parent_id
            ind_parents = [female_parent, male_parent]
            ind_time = max_gen - pedigree_df.iloc[x, 1]
            this_pedigree.add_individual(
                time = ind_time, 
                parents = ind_parents, 
                population = "H"
            )
        this_pedigree.add_individual(
            time = 0, 
            parents = [0, 0], 
            population = "S1")
        this_pedigree.add_individual(
            time = 0, 
            parents = [1, 1], 
            population = "S2")
            
        pedigree_ts = this_pedigree.finalise(sequence_length = self.seq_length)
        self.pedigree_ts = pedigree_ts
        
        return("pedigree successfully created")
        
    def make_sim(self):
        
        self.make_pedigree()
        sim_ts = msprime.sim_ancestry(
            initial_state = self.pedigree_ts,
            model = 'fixed_pedigree',
            ploidy = 2,
            recombination_rate = self.recomb_rate
        )
        self.sim_ts = sim_ts
        
        return("simulation successfully run")    
    
    def make_root_list(self):
        
        root_set = set()
        for tree_index in range(self.sim_ts.num_trees):
            tree = self.sim_ts.at_index(tree_index)
            roots = tree.roots
            root_set.update(roots)
        root_list = list(root_set)
        self.root_list = root_list
        
        return(root_list)
    
    def make_leaf_list(self):
        """xx"""
        tree = self.sim_ts.first()
        leaf_list = list(tree.leaves())
        
        return(leaf_list)
    
    def find_root(self, x, locus):
        """xx"""
        tree = self.sim_ts.at(locus)
        while tree.parent(x) != tskit.NULL:
            x = tree.parent(x)
            
        return(x)
    
    def make_squashed_edge(self):
        """ Makes the squashed edge table, a compressed tskit edge table which
            links each leaf to its population of origin.
            method:
                nodes (rows) whose children are not leaves (in the last 
                    generation) are removed
                nodes whose pants are not roots have their parents replaced by 
                    appropriate roots found using the find_root function
                parent is replaced with the provenance of the root (1 or 2)
                table is squashed to merge all adjacent tracts with the same 
                    provenance and order the table, then returned
        """
        edge_table = self.sim_ts.tables.edges
        leaf_list = self.make_leaf_list()
        root_list = self.make_root_list()
        fake_col = tskit.TableCollection()
        squashed_edge = fake_col.edges
        
        for row in edge_table:
            if row.child in leaf_list and self.sim_ts.tables.nodes[row.child].population == 0:
                squashed_edge.append(row)
        for n in range(len(squashed_edge)):
            row = squashed_edge[n] 
            if row.parent not in root_list:     
                locus = row.left 
                root = self.find_root(row.child, locus)
                squashed_edge[n] = squashed_edge[n].replace(parent = root)
        for n in range(len(squashed_edge)):
            prov = self.sim_ts.tables.nodes[squashed_edge[n].parent].population
            squashed_edge[n] = squashed_edge[n].replace(parent = prov)
        squashed_edge.squash()
        self.squashed_edge = squashed_edge
        
        return(squashed_edge)
    
    def make_genome_df(self):
        """ Runs the make_squashed_edge function to create a squashed_edge
            table, then transforms its columns into lists and integrates them
            into a dataframe under the columns
            Id
            Chr
            Start(bp)
            End(bp)
            Start(cM)
            End(cM)
            ans returns the dataframe.
        """
        squashed_edge = self.make_squashed_edge()
        
        Id = []
        Chr = []
        Start_bp = []
        End_bp = []
        Provs = []
        Start_cM = []
        End_cM = []
    
        for row in squashed_edge:
            Id.append(row.child // 2)
            Chr.append(row.child % 2)
            Start_bp.append(int(row.left))
            End_bp.append(int(row.right))
            Provs.append(row.parent)
            Start_cM.append(row.left * self.recomb_rate)
            End_cM.append(row.right * self.recomb_rate)
        
        table = {"Id": Id,
                 "Chr": Chr,
                 "Start(bp)": Start_bp,
                 "End(bp)": End_bp,
                 "Prov": Provs,
                 "Start(cM)": Start_cM,
                 "End(cM)": End_cM}
        
        genome_df = pd.DataFrame(table)
        genome_df = genome_df.sort_values(by = ["Id", "Chr", "Start(bp)"])
        genome_df = genome_df.reset_index(drop = True)
        
        return(genome_df)
    
    def write_genome_df(self, output_filename):
        """runs make_genome_df and writes the genome_df to file"""
        genome_df = self.make_genome_df()
        genome_string = genome_df.to_string(index = False)
        output_file = open(output_filename, 'w')
        output_file.write(genome_string)
        output_file.close()
        
        return("file output successful")
  
    def make_tree_svg(self):
        """ displays an svg tree scaled very roughly based on 
            sequence length
        """
        node_labels = {node.id: f"{node.id}" for node in self.sim_ts.nodes()}
        #{node.individual} before node.id
        sim_svg = self.sim_ts.draw_svg(
            size = (200 * self.seq_length, 400),
            node_labels=node_labels,
            y_axis = True,
            y_ticks = [0, 5, 10, 15, 20, 25, 30])
        display(SVG(sim_svg))

    def make_web(self):
        """ a utility for visualizing lineages using networkx. credit 
            to the msprime manual at
            https://tskit.dev/msprime/docs/latest/pedigrees.html
        """
        G = nx.DiGraph()
        for ind in self.sim_ts.individuals():
            time = self.sim_ts.node(ind.nodes[0]).time
            G.add_node(ind.id, time=time)
            for p in ind.parents:
                if p != tskit.NULL:
                    G.add_edge(ind.id, p)
        pos = nx.multipartite_layout(G, subset_key="time", 
            align="horizontal")
        fig, ax = plt.subplots(figsize=(20,20))
        nx.draw_networkx(G, pos, with_labels=True)
