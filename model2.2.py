#a new model, where I've changed the way trials are done by making them into class instances. Now many can be run at the same time and compared directly!

import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import copy

class Trial:
    
    def __init__(self, description, runs, g, S1f_init, S1m_init, S2f_init, S2m_init, 
                S1f_const, S1m_const, S2f_const, S2m_const,
                c11_0, c22_0, c11, chh, c22):
        
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
        self.c1h = -(c11/2)
        self.c12 = -(c11/2)
        
        self.ch1 = -(chh/2)
        self.chh = chh
        self.ch2 = -(chh/2)
        
        self.c21 = -(c22/2)
        self.c2h = -(c22/2)
        self.c22 = c22
        
        self.N_init = self.S1f_init + self.S1m_init + self.S2f_init + self.S2m_init
        self.N_const = self.S1f_const + self.S1m_const + self.S2f_const + self.S2m_const
        
        
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

        return(gen_0)
    

    def constant_contribution(self, gen_g_minus1):
        '''appends constant contributions to gen_g_minus1'''

        for x in range(self.S1f_const):
            gen_g_minus1.append([0, 1, 1])

        for x in range(self.S1m_const):
            gen_g_minus1.append([1, 1, 1])

        for x in range(self.S2f_const):
            gen_g_minus1.append([0, 0, 2])

        for x in range(self.S2m_const):
            gen_g_minus1.append([1, 0, 2])

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
    

    def mate(self, gen_g_minus1, is_gen_1):
        '''creates a new generation'''

        #define probabilities- gen 1 is a special case with its own prob function
        if is_gen_1 == True:
            prob_list = self.find_probs_gen1(gen_g_minus1)
        else:
            prob_list = self.find_probs(gen_g_minus1)

        #partition individuals into lists based on sex and population, adding their original index to them
        S1f_inds, S1m_inds, Hf_inds, Hm_inds, S2f_inds, S2m_inds = [], [], [], [], [], []

        #this is to avoid appending indexes to gen_g_minus1
        waste_list = copy.deepcopy(gen_g_minus1)

        for x in range(len(waste_list)):
            ind = waste_list[x]
            #x is the index of the individual within the generation
            ind.append(x)
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

        return(gen_g, gen_g_parents)
    

    def tabulate_parents_gen_0(self, gen_0):
        '''tabulates parentage for gen 0 individuals. all parents in gen 0 are -1 or -2'''

        ancestries = []
        pops = []
        female_parents = []
        male_parents = []

        for index in range(len(gen_0)):

            ind = gen_0[index]
            pops.append(ind[2])

            if ind[2] == 1:
                ancestries.append(1)
                female_parents.append(-1)
                male_parents.append(-1)

            elif ind[2] == 2:
                ancestries.append(0)
                female_parents.append(-2)
                male_parents.append(-2)

        gen_0_parent_table = {'g' : [0] * len(gen_0),
                              'index' : range(len(gen_0)),
                              'ancestry' : ancestries,
                              'pop' : pops,
                              'female_parent' : female_parents,
                              'male_parent' : male_parents   }
        gen_0_parent_frame = pd.DataFrame(gen_0_parent_table)

        return(gen_0_parent_frame)
    

    def tabulate_parents(self, gen_g, gen_g_parents, x):
        '''tabulates parentage for gen g individuals. x is generation number'''

        ancestries = []
        pops = []
        female_parents = []
        male_parents = []

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
        else:

            for index in range(len(gen_g)):

                parents = gen_g_parents[index]
                female_parent = parents[0]
                male_parent = parents[1]

                if female_parent[2] == 1:
                    female_parents.append(-1)    
                elif female_parent[2] == 2:
                    female_parents.append(-2)  
                else:
                    female_parents.append(female_parent[3])

                if male_parent[2] == 1: 
                    male_parents.append(-1)
                elif male_parent[2] == 2:
                    male_parents.append(-2)
                else:
                    male_parents.append(male_parent[3])

                ind = gen_g[index]
                pops.append(ind[2])
                ancestries.append(ind[1])

        gen_g_parent_table = {'g' : [x] * len(gen_g),
                              'index' : range(len(gen_g)),
                              'ancestry' : ancestries,
                              'pop' : pops,
                              'female_parent' : female_parents,
                              'male_parent' : male_parents   }
        gen_g_parent_frame = pd.DataFrame(gen_g_parent_table)

        return(gen_g_parent_frame)
    

    def make_gens(self):
        '''produces a list of generations'''
        gens = []
        parent_frame_list = []

        gen_0 = self.new_population()
        gens.append(gen_0) 
        gen_g_minus1 = copy.deepcopy(gen_0)
        gen_0_parent_frame = self.tabulate_parents_gen_0(gen_0)
        parent_frame_list.append(gen_0_parent_frame)

        for x in range(1, self.g + 1):

            if x == 1:
                gen_g, gen_g_parents = self.mate(gen_g_minus1, True)
                gens.append(gen_g)

            else:
                del(gen_g_minus1[(self.N_init - self.N_const):])
                gen_g_minus1 = self.constant_contribution(gen_g_minus1)
                gen_g, gen_g_parents = self.mate(gen_g_minus1, False)
                gens.append(gen_g)

            gen_g_parent_frame = self.tabulate_parents(gen_g, gen_g_parents, x)
            parent_frame_list.append(gen_g_parent_frame)

            gen_g_minus1 = copy.deepcopy(gen_g)

        parent_frames = pd.concat(parent_frame_list)

        return(gens, parent_frames)
    

    def execute_trial(self):

        many_gens = []
        many_parent_frames_list = []

        for x in range(self.runs):
            gens, parent_frames = self.make_gens()
            many_gens.append(gens)
            many_parent_frames_list.append(parent_frames)
        
        self.many_gens = many_gens
        self.many_parent_frames = pd.concat(many_parent_frames_list)
        
        
    def var_plot(self):
        var = []
        testgens = self.many_gens[0]
        for x in range(len(testgens)):
            gen_x_ancestries = []
            for gens in self.many_gens:
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

        gens = self.many_gens[0]
        for x in range(0, self.g):
            ancestries = []
            for y in gens[x]:
                ancestries.append(y[1]) 

            variance = round(np.var(ancestries), 5)
            
            ax[x].hist(ancestries, bins = 100, range = (0, 1), density = True, color = 'red')
            ax[x].set(title = ('gen = {}, variance = {}'.format(x, variance)))
            ax[x].set_xlim(-0.05, 1.05)
            ax[x].set_ylim(0, 100)


trial_1 = Trial(description = 'N = 2000, m = 800 random',
               runs=1, g=20, 
               S1f_init=500, S1m_init=500, S2f_init=500, S2m_init=500, 
               S1f_const=200, S1m_const=200, S2f_const=200, S2m_const=200,
               c11_0=0, c22_0=0, 
               c11=0, chh=0, c22=0)

trial_2 = Trial(description = 'N = 2000, m = 800 c11=c22=0.04',
               runs=1, g=20, 
               S1f_init=500, S1m_init=500, S2f_init=500, S2m_init=500, 
               S1f_const=200, S1m_const=200, S2f_const=200, S2m_const=200,
               c11_0=0, c22_0=0, 
               c11=0.04, chh=0, c22=0.04)

trial_3 = Trial(description = 'N = 2000, m = 800 c11=c22=-0.04',
               runs=1, g=20, 
               S1f_init=500, S1m_init=500, S2f_init=500, S2m_init=500, 
               S1f_const=200, S1m_const=200, S2f_const=200, S2m_const=200,
               c11_0=0, c22_0=0, 
               c11=-0.04, chh=0, c22=-0.04)

trial_1.execute_trial()
trial_2.execute_trial()
trial_3.execute_trial()

varfig = plt.figure()
varx = varfig.add_subplot(111)
varx.set(xlabel = 'generations', ylabel = 'variance of S1 ancestry',
              title = 'Variance in Ancestry')

trial_1.var_plot()
trial_2.var_plot()
trial_3.var_plot()

varx.legend([trial_1.description, trial_2.description, trial_3.description])
