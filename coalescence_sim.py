import tskit
import msprime
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import SVG
import networkx as nx
import os


class CoalescenceSim:
   
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
        for f in os.listdir(output_dir):
            os.remove(os.path.join(output_dir, f))

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