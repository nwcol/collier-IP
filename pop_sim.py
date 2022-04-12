import random
import pandas as pd
import copy


class PopSim:
    
    def __init__(
        self,
        g, 
        seed = None,
        runs = 1,  
        S1f_init = 0, 
        S1m_init = 0, 
        S2f_init = 0, 
        S2m_init = 0, 
        N_init = None,
        S1_init = None,
        S2_init = None,
        S1f_const = 0, 
        S1m_const = 0, 
        S2f_const = 0, 
        S2m_const = 0,
        N_const =None,
        c11_0 = 0,  
        c11 = 0, 
        c1h = 0, 
        c12 = 0,
        ch1 = 0, 
        chh = 0, 
        ch2 = 0,
        c21 = 0, 
        c2h = 0, 
        c22 = 0):
        """ Initialize the model.
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
        if seed == None:
            self.seed = self.random_seed()
        else:
            seld.seed = seed
        random.seed(self.seed)
        self.runs = runs
        self.g = g
        
        if N_init == None and S1_init == None:
            self.S1f_init = S1f_init
            self.S1m_init = S1m_init
            self.S2f_init = S2f_init
            self.S2m_init = S2m_init
        elif S1_init != None:
            self.S1f_init = S1_init // 2
            self.S1m_init = S1_init // 2
            self.S2f_init = S2_init // 2
            self.S2m_init = S2_init // 2
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
        self.c22_0 = c11_0
        self.c21_0 = -c11_0
        
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
    
    def random_seed(self):
        num = range(0, 10)
        n = 8
        num_list = random.sample(num, n)
        if len(num_list) < 8:
            num_list.append(random.sample(num))
        seed = [str(x) for x in num_list]
        seed = "".join(seed)
        seed = int(seed)
        return(seed)
        
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
        """ appends constant migrations to gen_g_minus1"""
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
        """ tabulates parentage for gen 0 individuals. all parents in 
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
        """ produces a list of generations"""        
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
        """ Integrates functions to perform a series of trials"""
        gens_list = []
        pedigree_df_list = []
        for x in range(self.runs):
            gens, pedigree_df = self.make_gens()
            gens_list.append(gens)
            pedigree_df_list.append(pedigree_df)
        self.gens_list = gens_list
        self.pedigree_df_list = pedigree_df_list 
    
    def write_pedigree_df(self, output_filename, df_index = 0):
        selected_df = self.pedigree_df_list[df_index]
        selected_df.to_csv(output_filename) 
        return("file output complete")

    def sample_pedigree_df(self, sample_size, df_index = 0):
        """ samples a specificed number of individuals from the last 
            generation and trims all individuals which aren't 
            ancestral to them out"""
        pedigree_df = self.pedigree_df_list[df_index]
        #pedigree_df = pedigree_df.drop(columns=["Unnamed: 0"])
        sample_df = pedigree_df.loc[pedigree_df['g'] == self.g].sample(n = sample_size)
        thisgen_df = sample_df
        for n in range(self.g):
            fem = thisgen_df["female_parent"].tolist()
            male = thisgen_df["male_parent"].tolist()
            parents = list(set(fem + male))
            if -1 in parents:
                parents.remove(-1)
            if -2 in parents:
                parents.remove(-2)
            lastgen_df = pedigree_df.iloc[parents]
            con = sample_df, lastgen_df
            sample_df = pd.concat(con)
            thisgen_df = lastgen_df
        sample_df = sample_df.sort_values(by = ["g", "id"])

        sample_df["g"].apply(lambda x: abs(x - self.g))
        sample_df["g"] = sample_df["g"].astype("string")
        for row in range(len(sample_df)):
            replaced_id = sample_df.iloc[row, 1]
            sample_df = sample_df.replace(to_replace = replaced_id, value = str(row))
        sample_df["g"] = pd.to_numeric(sample_df["g"])
        sample_df = sample_df.reset_index(drop = True)
        return(sample_df)
      
    def write_sample_df(self, sample_size, output_filename, df_index = 0):
        """ write sample_df to file as .csv"""
        sample_df = self.sample_pedigree_df(sample_size)
        sample_df.to_csv(output_filename) 
        return("file output complete")