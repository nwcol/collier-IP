import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import copy

class Model:
    
    
    def __init__(self, description, runs, g, S1f_init, S1m_init, S2f_init, S2m_init, 
                S1f_const, S1m_const, S2f_const, S2m_const,
                c11_0, c22_0, 
                c11, c1h, c12,
                ch1, chh, ch2,
                c21, c2h, c22):
        
        self.description = description
    
        self.runs = runs
        self.g = g
        
        self.S1f_init = S1f_init
        self.S1m_init = S1m_init
        self.S2f_init = S2m_init
        self.S2m_init = S2f_init
        
        self.S1f_const = S1f_const
        self.S1m_const = S1m_const
        self.S2f_const = S2f_const
        self.S2m_const = S2m_const
        
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
        
        self.N_init = self.S1f_init + self.S1m_init + self.S2f_init + self.S2m_init
        self.N_const = self.S1f_const + self.S1m_const + self.S2f_const + self.S2m_const
        
        self.N_f_init = self.S1f_init + self.S2f_init
        self.N_m_init = self.S1m_init + self.S2m_init
        
        self.id_shift = 0
        
        
    def new_population(self):
        '''generates an initial generation 0 based on initial parameters'''

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
        '''appends constant migrations to gen_g_minus1'''

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
        '''returns probabilities of matings for gen 1'''

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
        '''returns probabilities of matings based on composition of gen g-1'''

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

        #interim way to remove negative probabilities
        prob_list = [P11, P1h, P12, Ph1, Phh, Ph2, P21, P2h, P22]
        #for x in range(9):
         #   if prob_list[x] < 0:
         #       prob_list[x] = 0          

        return(prob_list)
    

    def mate(self, gen_g_minus1, is_gen_1, g):
        '''creates a new generation'''

        if is_gen_1 == True:
            prob_list = self.find_probs_gen1(gen_g_minus1)
        else:
            prob_list = self.find_probs(gen_g_minus1)

        S1f_inds, S1m_inds, Hf_inds, Hm_inds, S2f_inds, S2m_inds = [], [], [], [], [], []
        
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

        #the actual mating element
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
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '1h':
                female_parent = random.choice(S1f_inds)
                male_parent = random.choice(Hm_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '12':
                female_parent = random.choice(S1f_inds)
                male_parent = random.choice(S2m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == 'h1':
                female_parent = random.choice(Hf_inds)
                male_parent = random.choice(S1m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == 'hh':
                female_parent = random.choice(Hf_inds)
                male_parent = random.choice(Hm_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == 'h2':
                female_parent = random.choice(Hf_inds)
                male_parent = random.choice(S2m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '21':
                female_parent = random.choice(S2f_inds)
                male_parent = random.choice(S1m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '2h':
                female_parent = random.choice(S2f_inds)
                male_parent = random.choice(Hm_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring)

            elif mating == '22':
                female_parent = random.choice(S2f_inds)
                male_parent = random.choice(S2m_inds)
                gen_g_parents.append([female_parent, male_parent])
                offspring = [random.randint(0, 1), (female_parent[1] + male_parent[1]) / 2, 0]
                gen_g.append(offspring) 
        
        for x in range(len(gen_g)):
            gen_g[x].append(self.id_shift + x)
        
        self.id_shift = self.id_shift + len(gen_g)

        return(gen_g, gen_g_parents)
    

    def tabulate_pedigree_gen_0(self, gen_0):
        '''tabulates parentage for gen 0 individuals. all parents in gen 0 are -1 or -2'''

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
        '''tabulates parentage for gen g individuals. x is generation number'''

        ids = []
        ancestries = []
        pops = []
        female_parents = []
        male_parents = []
        
        migrant_ids = []
        migrant_ancestries = []
        migrant_pops = []
        female_migrant_parents = []
        male_migrant_parents = []
        
        start_id = self.id_shift
        
        gen_g_length = len(gen_g)
        
        if x == 1:
            
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
                
            migrant_length = 0
        
        elif x > 1:

            for index in range(len(gen_g)):

                parents = gen_g_parents[index]
                female_parent = parents[0]
                male_parent = parents[1]

                ind = gen_g[index]
                pops.append(ind[2])
                ancestries.append(ind[1])
                ids.append(ind[3])

                if female_parent[2] != 0:

                    if female_parent[2] == 1:

                        migrant_ids.append(self.id_shift - gen_g_length)
                        migrant_ancestries.append(1)
                        migrant_pops.append(1)
                        female_migrant_parents.append(-1)
                        male_migrant_parents.append(-1)

                        self.id_shift = self.id_shift + 1

                        female_parents.append(migrant_ids[-1])

                    elif female_parent[2] == 2:

                        migrant_ids.append(self.id_shift - gen_g_length)
                        migrant_ancestries.append(0)
                        migrant_pops.append(2)
                        female_migrant_parents.append(-2)
                        male_migrant_parents.append(-2)

                        self.id_shift = self.id_shift + 1

                        female_parents.append(migrant_ids[-1])

                elif female_parent[2] == 0:

                    female_parents.append(female_parent[3])

                if male_parent[2] != 0:

                    if male_parent[2] == 1:

                        migrant_ids.append(self.id_shift - gen_g_length)
                        migrant_ancestries.append(1)
                        migrant_pops.append(1)
                        female_migrant_parents.append(-1)
                        male_migrant_parents.append(-1)

                        self.id_shift = self.id_shift + 1

                        male_parents.append(migrant_ids[-1])

                    elif male_parent[2] == 2:

                        migrant_ids.append(self.id_shift - gen_g_length)
                        migrant_ancestries.append(0)
                        migrant_pops.append(2)
                        female_migrant_parents.append(-2)
                        male_migrant_parents.append(-2)

                        self.id_shift = self.id_shift + 1

                        male_parents.append(migrant_ids[-1])

                elif male_parent[2] == 0:

                    male_parents.append(male_parent[3])


                migrant_length = len(migrant_ids)
        
        for z in range(len(ids)):

            ids[z] = ids[z] + migrant_length
            

        migrant_pedigree_table = {'g' : [x - 1] * migrant_length,
                               'id' : migrant_ids,
                               'ancestry' : migrant_ancestries,
                               'pop' : migrant_pops,
                               'female_parent' : female_migrant_parents,
                               'male_parent' : male_migrant_parents   }

        gen_g_pedigree_table = {'g' : [x] * len(gen_g),
                              'id' : ids,
                              'ancestry' : ancestries,
                              'pop' : pops,
                              'female_parent' : female_parents,
                              'male_parent' : male_parents   }
        
        gen_g_migrant_df = pd.DataFrame(migrant_pedigree_table)
        gen_g_pedigree_df_base = pd.DataFrame(gen_g_pedigree_table)
        
        gen_g_pedigree_df = pd.concat([gen_g_migrant_df, gen_g_pedigree_df_base])
        
        gen_g_pedigree_df = gen_g_pedigree_df.astype({
            'g': 'int',
            'id': 'int',
            'pop': 'int',
            'female_parent': 'int',
            'male_parent': 'int'
        })

        return(gen_g_pedigree_df)
    

    def make_gens(self):
        '''produces a list of generations'''
        
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

            gen_g_pedigree_df = self.tabulate_pedigree(gen_g, gen_g_parents, x)
            internal_pedigree_list.append(gen_g_pedigree_df)

            gen_g_minus1 = copy.deepcopy(gen_g)

        pedigree_df = pd.concat(internal_pedigree_list)

        return(gens, pedigree_df)
    

    def execute_model(self):

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
            
            ax[x].hist(ancestries, bins = 100, range = (0, 1), density = True, color = 'red')
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
        hf = (self.N_f_init - self.S1f_const - self.S2f_const) / self.N_f_init
        hm = (self.N_m_init - self.S1m_const - self.S2m_const) / self.N_m_init
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
                E_H_list.append((s1f_0 * s1m_0) + c11_0 + (((s1f_0 * s2m_0) + (s2f_0 * s1m_0) + c12_0 + c21_0) / 2))
                #from eq 20
                predicted_var.append(((s1f_0 * (1 - s1f_0)) + (s1m_0 * (1 - s1m_0)) + (2 * c11_0)) / 4)

            elif g > 1:
                #from eq 17
                E_H_list.append(((s1f * s1m) + c11) +
                               (((s1f * s2m) + (s2f * s1m) + c12 + c21 + c1h + ch1) / 2) +
                               ((((s1f * hm) + (hf * s1m) + c1h + ch1) / 2) * (1 + E_H_list[g-1])) +
                               ((((hf * hm) + chh) / 2) * (2 * E_H_list[g-1])) +
                               ((((s2f * hm) + (hf * s2m) + ch2 + c2h) / 2) * E_H_list[g-1])
                               )                     
                
                #from eq 21
                predicted_var.append((((s1f * (1 - s1f)) + (s1m * (1 - s1m)) + (2 * c11)) / 4) +
                                    (((c1h + ch1 - (s1f * hf) - (s1m * hm)) / 2) * E_H_list[g-1]) +
                                    ((((hf * (1 - hf)) + (hm * (1 - hm)) + (2 * chh)) / 4) * (E_H_list[g-1] ** 2)) +
                                    (((hf + hm) / 4) * predicted_var[g-1])
                                    )
        
        x = range(len(var))
        predicted_varfig = plt.figure()
        var_predx = predicted_varfig.add_subplot(111)
        var_predx.plot(x, var, color = ('red'))
        var_predx.plot(predicted_var, color = ('black'))
        var_predx.set_ylim(0, .25)
        var_predx.set_xlim(0, len(x))
        plt.xticks(np.arange(min(x), max(x) + 1, 5))
        
        var_predx.set(xlabel = 'generations', ylabel = 'variance of S1 ancestry', title = 'Variance in Ancestry')
        var_predx.legend([self.description, 'predicted variance'])
        
            
    
    def write_pedigree_df(self, df_index, output_filename):
        
        selected_df = self.pedigree_df_list[df_index]
        selected_df.to_csv(output_filename) 
        return("file output successful")
        
    
    def write_pedigree_df_list(self, output_filename = None):
        
        if output_filename == None:
            
            output_filename = "model df_list N" + str(self.N_init) + " m" + str(self.N_const) + ".csv"
            output_filename = "modeldata/" + output_filename     
            
        output_file = open(output_filename, 'w')

        for pedigree_df in self.pedigree_df_list:
            pedigree_df.to_csv(output_file)
            
        output_file.close
        return("file output successful")
            
        

model_1 = Model(description = 'c11 = c22 = 0', 
               runs=1, g=10, 
               S1f_init=50, S1m_init=50, S2f_init=50, S2m_init=50, 
               S1f_const=0, S1m_const=0, S2f_const=0, S2m_const=0,
               c11_0=0, c22_0=0,
               c11=0, c1h=0, c12=0,
               ch1=0, chh=0, ch2=0,
               c21=0, c2h=0, c22=0)

model_1.execute_model()
model_1.write_pedigree_df(df_index = 0, output_filename = "written_parent_frame.csv")
