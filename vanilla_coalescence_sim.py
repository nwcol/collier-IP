import numpy as np
import tskit
import msprime
import pandas as pd
import os
import sys
tractspath = "C:/Genetics work/tracts-python3"
sys.path.append(tractspath)
import tracts

class VanillaCoalescenceSim:
    
    def __init__(self, admix_time, pop_size, sample_size, cM_length, seq_length, s1_init = 0.5, s2_init = 0.5):
        self.admix_time = admix_time
        self.s1_init = s1_init
        self.s2_init = s2_init
        self.pop_size = pop_size
        self.sample_size = sample_size
        self.cM_length = cM_length
        self.seq_length = seq_length
        self.recomb_rate = cM_length * 0.01 / seq_length

    def execute_sim(self, output_dir):
        vanilla_ts = self.coalescent_sim()
        genome_df = self.make_genome_df(vanilla_ts)
        self.write_genomes(genome_df, output_dir)
        
    def coalescent_sim(self):
        demog = msprime.Demography()
        demog.add_population(name = "H", initial_size = self.pop_size)
        demog.add_population(name = "S1", initial_size = 10)
        demog.add_population(name = "S2", initial_size = 10)
        demog.add_population(name = "ANC", initial_size = 100)
        demog.add_admixture(
            time = self.admix_time, 
            derived = "H", 
            ancestral = ["S1", "S2"], 
            proportions = [self.s1_init, self.s2_init])
        demog.add_population_split(time = self.admix_time + 100, derived=["S1", "S2"], ancestral="ANC")
        sample_set = {"H": self.sample_size, "S1": 1, "S2": 1}
        
        vanilla_ts = msprime.sim_ancestry(
            model = "dtwf",
            demography = demog,
            samples = sample_set,
            ploidy = 2,
            recombination_rate = self.recomb_rate,
            sequence_length = self.seq_length)
        self.vanilla_ts = vanilla_ts
        return(vanilla_ts)
    
    def make_squashed_edge(self, vanilla_ts):
        leaf_list = self.make_leaf_list(vanilla_ts)
        fake_edge = self.remove_anc(vanilla_ts)
        squashed_col = tskit.TableCollection()
        squashed_edge = squashed_col.edges
        
        for row in vanilla_ts.tables.edges:
            if row.child in leaf_list and vanilla_ts.tables.nodes[row.child].population == 0:
                squashed_edge.append(row)            
        for n in range(len(squashed_edge)):
            row = squashed_edge[n]     
            locus = row.left 
            root = self.find_root(row.child, locus)
            squashed_edge[n] = squashed_edge[n].replace(parent = root)
        for n in range(len(squashed_edge)):
            prov = vanilla_ts.tables.nodes[squashed_edge[n].parent].population
            squashed_edge[n] = squashed_edge[n].replace(parent = prov)
        squashed_edge.squash()
        return(squashed_edge)
    
    def remove_anc(self, vanilla_ts):
        #make a new edge table without members of population 3
        anc_roots = []
        for n in range(len(vanilla_ts.tables.nodes)):
            if vanilla_ts.tables.nodes[n].population == 3:
                anc_roots.append(n)
        fake_col = tskit.TableCollection()
        fake_edge = fake_col.edges
        for row in vanilla_ts.tables.edges:
            if row.parent not in anc_roots:
                fake_edge.append(row)
        return(fake_edge)
    
    def make_leaf_list(self, vanilla_ts):
        """xx"""
        tree = self.vanilla_ts.first()
        leaf_list = list(tree.leaves())
        return(leaf_list)
    
    def find_root(self, x, locus):
        """xx"""
        tree = self.vanilla_ts.at(locus)
        while self.vanilla_ts.tables.nodes[tree.parent(x)].population != 3:
            x = tree.parent(x)
        return(x)

    def make_genome_df(self, vanilla_ts):
        squashed_edge = self.make_squashed_edge(vanilla_ts)
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
    
    def write_genomes(self, genome_df, output_dir):
        """ runs make_genome_df and writes each chromosome to its own file"""
        #delete any folders lying around in the output folder
        #kill_files = os.listdir(output_dir)
        #for f in kill_files:
        #   os.remove(f)
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