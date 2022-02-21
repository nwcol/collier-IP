import numpy as np
import matplotlib.pyplot as plt
import random
import pandas as pd
import copy

class Trial:
    
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
            
        return(gen_0)
    

    def constant_contribution(self, gen_g_minus1):
        '''appends constant contributions to gen_g_minus1'''

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
                
        start_id = self.N_init * g
        for x in range(len(gen_g)):
            gen_g[x].append(start_id + x)

        return(gen_g, gen_g_parents)
    

    def tabulate_parents_gen_0(self, gen_0):
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

        gen_0_parent_table = {'g' : [0] * len(gen_0),
                              'id' : range(len(gen_0)),
                              'ancestry' : ancestries,
                              'pop' : pops,
                              'female_parent' : female_parents,
                              'male_parent' : male_parents   }
        gen_0_parent_frame = pd.DataFrame(gen_0_parent_table)

        return(gen_0_parent_frame)
    

    def tabulate_parents(self, gen_g, gen_g_parents, x):
        '''tabulates parentage for gen g individuals. x is generation number'''

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

        gen_g_parent_table = {'g' : [x] * len(gen_g),
                              'id' : ids,
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
                gen_g, gen_g_parents = self.mate(gen_g_minus1, True, x)
                gens.append(gen_g)

            else:
                del(gen_g_minus1[(self.N_init - self.N_const):])
                gen_g_minus1 = self.constant_contribution(gen_g_minus1)
                gen_g, gen_g_parents = self.mate(gen_g_minus1, False, x)
                gens.append(gen_g)

            gen_g_parent_frame = self.tabulate_parents(gen_g, gen_g_parents, x)
            parent_frame_list.append(gen_g_parent_frame)

            gen_g_minus1 = copy.deepcopy(gen_g)

        parent_frames = pd.concat(parent_frame_list)

        return(gens, parent_frames)
    

    def execute_trial(self):

        many_gens = []
        parent_frame_list = []

        for x in range(self.runs):
            gens, parent_frames = self.make_gens()
            many_gens.append(gens)
            parent_frame_list.append(parent_frames)
        
        self.many_gens = many_gens
        self.parent_frame_list = parent_frame_list
        
        
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
            
        
    def compare_predicted_var(self):
        
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
        
            
    
    def write_parent_frame(self, frame_index):
        selected_frame = self.parent_frame_list[frame_index]
        selected_frame.to_csv("written_parent_frame.csv") 
        
        

trial_1 = Trial(description = 'c11 = c22 = 0', 
               runs=100, g=20, 
               S1f_init=500, S1m_init=500, S2f_init=500, S2m_init=500, 
               S1f_const=200, S1m_const=200, S2f_const=200, S2m_const=200,
               c11_0=0, c22_0=0,
               c11=0.02, c1h=-0.04, c12=0.02,
               ch1=-0.04, chh=0.02, ch2=0.02,
               c21=0.02, c2h=0.02, c22=-0.04))

trial_2 = Trial(description = 'c11 = c22 = 0.02',
               runs=100, g=20, 
               S1f_init=500, S1m_init=500, S2f_init=500, S2m_init=500, 
               S1f_const=200, S1m_const=200, S2f_const=200, S2m_const=200,
               c11_0=0, c22_0=0, 
               c11=0.02, c1h=0, c12=-0.02,
               ch1=0, chh=0, ch2=0,
               c21=-0.02, c2h=-0, c22=0.02)

trial_3 = Trial(description = 'c11 = c22 = -0.04',
               runs=100, g=20, 
               S1f_init=500, S1m_init=500, S2f_init=500, S2m_init=500, 
               S1f_const=200, S1m_const=200, S2f_const=200, S2m_const=200,
               c11_0=0, c22_0=0, 
               c11=-0.04, c1h=0, c12=0.04,
               ch1=0, chh=0, ch2=0,
               c21=0.04, c2h=0, c22=-0.04)

trial_4 = Trial(description = 'c11 = c22 = -0.02', 
               runs=100, g=20, 
               S1f_init=500, S1m_init=500, S2f_init=500, S2m_init=500, 
               S1f_const=200, S1m_const=200, S2f_const=200, S2m_const=200,
               c11_0=0, c22_0=0,
               c11=-0.02, c1h=0, c12=0.02,
               ch1=0, chh=0, ch2=0,
               c21=0.02, c2h=0, c22=-0.02)
                
trial_5 = Trial(description = 'c11 = c22 = 0.04', 
               runs=100, g=20, 
               S1f_init=500, S1m_init=500, S2f_init=500, S2m_init=500, 
               S1f_const=200, S1m_const=200, S2f_const=200, S2m_const=200,
               c11_0=0, c22_0=0,
               c11=0.04, c1h=0, c12=-0.04,
               ch1=0, chh=0, ch2=0,
               c21=-0.04, c2h=0, c22=0.04)

trial_1.execute_trial()
trial_2.execute_trial()
trial_3.execute_trial()
trial_4.execute_trial()
trial_5.execute_trial()

#trial_1.write_parent_frame(0)

varfig = plt.figure()
varx = varfig.add_subplot(111)
varx.set(xlabel = 'generations', ylabel = 'variance of S1 ancestry',
              title = 'Variance in Ancestry')
trial_1.var_plot()
trial_2.var_plot()
trial_3.var_plot()
trial_4.var_plot()
trial_5.var_plot()
varx.legend([trial_1.description, trial_2.description, trial_3.description, trial_4.description, trial_5.description])


trial_1.compare_predicted_var()
trial_2.compare_predicted_var()
trial_3.compare_predicted_var()
trial_4.compare_predicted_var()
trial_5.compare_predicted_var()
