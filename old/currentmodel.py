import numpy as np
import tskit
import msprime
import matplotlib.pyplot as plt
import random
import pandas as pd
import copy
from IPython.display import SVG
import networkx as nx
import os
import sys
import tkinter as Tk
import pylab

tractspath = "C:/Genetics work/tracts-python3"
sys.path.append(tractspath)
import tracts

class PopModel:
    
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
        
    def set_c(self):
        """calculates all c values using the input cs
            under construction!
        """
        pass
    
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
            fem = list(set(thisgen_df["female_parent"].tolist()))
            male = list(set(thisgen_df["male_parent"].tolist()))
            parents = fem + male
            lastgen_df = pedigree_df.iloc[parents]
            con = sample_df, lastgen_df
            sample_df = pd.concat(con)
            thisgen_df = lastgen_df
        sample_df = sample_df.sort_values(by = ["g", "id"])

        sample_df["g"].apply(lambda x: abs(x - self.g))
        sample_df["g"] = sample_df["g"].astype("string")
        for row in range(len(sample_df)):
            replaced_id = sample_df.iloc[row, 1]
            sample_df =sample_df.replace(to_replace = replaced_id, value = str(row))
        sample_df["g"] = pd.to_numeric(sample_df["g"])
        sample_df = sample_df.reset_index(drop = True)
        return(sample_df)
      
    def write_sample_df(self, sample_size, output_filename, df_index = 0):
        """ write sample_df to file as .csv"""
        sample_df = self.sample_pedigree_df(sample_size, df_index)
        sample_df.to_csv(output_filename) 
        return("file output complete")

    
class CoalescentSim:
   
    def __init__(self, seed, input_filename, seq_length, cM_length):
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
        self.recomb_rate = cM_length * 0.01 / seq_length
        self.seed = seed
        
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
            population = 2)
        this_pedigree.add_individual(
            time = max_gen + 1, 
            parents = [-1, -1], 
            population = 1)
        
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
            population = 2)
        this_pedigree.add_individual(
            time = 0, 
            parents = [1, 1], 
            population = 1)
            
        pedigree_ts = this_pedigree.finalise(sequence_length = self.seq_length)
        self.pedigree_ts = pedigree_ts
        
        return("pedigree successfully created")
        
    def make_sim(self):
        
        self.make_pedigree()
        sim_ts = msprime.sim_ancestry(
            initial_state = self.pedigree_ts,
            model = 'fixed_pedigree',
            ploidy = 2,
            recombination_rate = self.recomb_rate,
            random_seed = self.seed
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
            and returns the dataframe.
        """
        squashed_edge = self.make_squashed_edge()
        bp_to_M = 100 * self.recomb_rate
        
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
            Start_cM.append(row.left * bp_to_M)
            End_cM.append(row.right * bp_to_M)
        
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
        """ runs make_genome_df and writes the genome_df to file"""
        genome_df = self.make_genome_df()
        genome_string = genome_df.to_string(index = False)
        output_file = open(output_filename, 'w')
        output_file.write(genome_string)
        output_file.close()
        
        return("file output complete")
    
    def make_ind_list(self, genome_df):
        """ pulls the ids of each individual in a genome_df out of the genome_df
            into a list.
        """
        ind_list = list(set(genome_df['Id'].to_list()))
        return(ind_list)
    
    def partition_genome_df(self, genome_df, ind_list):
        """ partitions the genome_df into a list of chromosomes"""
        chrom_list = []
        for ind in ind_list:
            ind_df = genome_df.groupby(['Id']).get_group(ind)
            for chrom in range(0, 2):
                chrom_df = ind_df.groupby(['Chr']).get_group(chrom)  
                chrom_list.append(chrom_df)
        for x in range(len(chrom_list)):
            chrom_list[x] = chrom_list[x].drop(columns=['Id'])
            (chrom_list[x])['Chr'] = 'Chr1'
        return(chrom_list)
    
    def write_genomes(self, output_dir):
        """ runs make_genome_df and writes each chromosome to its own file"""
        #delete any folders lying around in the output folder
        #kill_files = os.listdir(output_dir)
        #for f in kill_files:
        #   os.remove(f)
        genome_df = self.make_genome_df()
        ind_list = self.make_ind_list(genome_df)
        chrom_list = self.partition_genome_df(genome_df, ind_list)
        file_names = []
        for x in range(len(ind_list)):
            file_names.append(output_dir + '/' + str(x) + '_A.bed')
            file_names.append(output_dir + '/' + str(x) + '_B.bed')
        for x in range(len(file_names)):
            chrom_string = chrom_list[x].to_string(index = False)
            output_file = open(file_names[x], 'w')
            output_file.write(chrom_string)
            output_file.close()
        return("file output complete")
  
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
        

class TwoParamInference:
    
    def __init__(self, input_dir):
        self.input_dir = input_dir
        
    def make_inference(self):
        pop, bins, data, labels = self.init_pop(self.input_dir)
        mig, xopt = self.optimize(pop, bins, data, labels)
        
        (t_start, s1_init) = xopt
        gen = t_start * 100
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
        
        mig = self.two_param_model(params) #get the migration matrix
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
