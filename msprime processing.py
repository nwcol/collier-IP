import numpy as np
import pandas as pd
import tskit
import msprime
from IPython.display import SVG
import json
import matplotlib.pyplot as plt

class Sim:
    
    
    def __init__(
                self,
                input_filename,
                seq_length,
                recomb_rate
                ):
    
        self.input_filename = input_filename
        self.seq_length = seq_length
        self.recomb_rate = recomb_rate
        

    def execute_sim(self):

        pedigree_df = pd.read_csv(self.input_filename)
        max_gen = pedigree_df.iloc[len(pedigree_df) - 1, 1]
        
        this_pedigree = msprime.PedigreeBuilder()

        pedigree_schema = tskit.MetadataSchema({
            "codec": "json",
            "root": {
                "simpleTypes": {
                  "enum": ["array", "boolean", "integer", "null", "number", "object", "string",]
                               }
                    },
            "properties": {
                "$comment": {"type": "string"}
                          }
        })

        this_pedigree.tables.individuals.metadata_schema = pedigree_schema

        for x in range(len(pedigree_df)):

            female_parent_id = pedigree_df.iloc[x, 5]
            male_parent_id = pedigree_df.iloc[x, 6]

            pop = str(pedigree_df.iloc[x, 4])

            if female_parent_id == -1:

                female_parent = -1

            elif female_parent_id == -2:

                female_parent = -1

            else:
                
                female_parent = female_parent_id

            if male_parent_id == -1:

                male_parent = -1

            elif male_parent_id == -2:

                male_parent = -1

            else:

                male_parent = male_parent_id
   
            ind_parents = [female_parent, male_parent]
            ind_time = max_gen - pedigree_df.iloc[x, 1]

            this_pedigree.add_individual(
                time=ind_time, 
                parents=ind_parents, 
                population=0,
                metadata={"$comment": pop}
            )

        pedigree_ts = this_pedigree.finalise(sequence_length = self.seq_length)

        sim_ts = msprime.sim_ancestry(
            initial_state = pedigree_ts,
            model = 'fixed_pedigree',
            ploidy = 2,
            recombination_rate = self.recomb_rate
        )
        
        self.sim_ts = sim_ts
                          
    
    def write_sim(self, 
                  output_filename):
                          
        output_file = open(output_filename, "w")
        self.sim_ts.dump(output_file)
        
    
    def show_tree(self):
    
        sim_svg = self.sim_ts.draw_svg(
            size = (200 * self.seq_length, 400),
            y_axis = True,
            y_ticks = [0, 5, 10, 15, 20, 25, 30]
        )
        display(SVG(sim_svg))
            
            
    def extract_genome_df_list(self):
        
        #select all individuals in time 0 and make list of their ids
        sample_inds = []
        sample_nodes = []

        for x in range(len(self.sim_ts.tables.nodes)):

            if x % 2 == 0:

                pass

            elif x % 2 == 1:

                selected_row = self.sim_ts.tables.nodes[x]

                if selected_row.time != 0:

                    pass

                elif selected_row.time == 0:

                    ind_id = selected_row.individual
                    sample_inds.append(ind_id)
                    nodes = [x - 1, x]
                    sample_nodes.append(nodes)

        genome_df_list = []

        #read edges table for gen 0
        for ind_index in range(len(sample_inds)):

            ind_chrom_list = []
            ind_start_list = []
            ind_stop_list = []
            ind_parent_node_list = []

            nodes = sample_nodes[ind_index]

            for node_index in range(2):

                node = nodes[node_index]

                for row in self.sim_ts.tables.edges:

                    if row.child == node:

                        ind_chrom_list.append(node_index)
                        ind_parent_node_list.append(row.parent)
                        ind_start_list.append(row.left)
                        ind_stop_list.append(row.right)

            #iterate through the parents list and find genome sources for each parent
            ind_source_list = []

            for parent_node in ind_parent_node_list:

                parent_pop = 0

                while parent_pop == 0:

                    parent_id = parent_node // 2
                    parent_pop = int(self.sim_ts.tables.individuals[parent_id].metadata['$comment'])

                    if parent_pop != 0:

                        ind_source_list.append(parent_pop)

                    else:

                        child_node = parent_node

                        for row in self.sim_ts.tables.edges:

                            if row.child == child_node:

                                parent_node = row.parent

            ind_genome = {"id": [sample_inds[ind_index]] * len(ind_chrom_list),
                          "chromosome": ind_chrom_list,
                          "starts": ind_start_list,
                          "stops": ind_stop_list,
                          "source": ind_source_list        
                        }

            ind_genome_df = pd.DataFrame(ind_genome)
            ind_genome_df = ind_genome_df.sort_values(by = ["chromosome", "starts"])
            genome_df_list.append(ind_genome_df)

        self.genome_df_list = genome_df_list

    
    def vis_genome_df_list(self):

        figure, ax = plt.subplots(figsize = (10, 0.25 * len(self.genome_df_list)))
        ax.set_ylim(0, len(self.genome_df_list))
        ax.set_xlim(0, self.seq_length)
        plt.grid(visible = True, axis = 'y', color = 'black', )

        for ind_genome in range(len(self.genome_df_list)):

            genome_table = self.genome_df_list[ind_genome]

            majtick = np.arange(0, len(self.genome_df_list), 1)
            mintick = np.arange(0, len(self.genome_df_list), 0.5)
            ax.set_yticks(majtick)
            ax.set_yticks(mintick, minor=True)

            labels = ['0', '1']

            chrom0_xranges = []
            chrom1_xranges = []
            chrom0_source = []
            chrom1_source = []

            for x in range(len(genome_table)):

                chrom_xrange = (genome_table.iloc[x, 2], 
                                genome_table.iloc[x, 3] - genome_table.iloc[x, 2])
                source = genome_table.iloc[x, 4]

                if genome_table.iloc[x, 1] == 0:

                    chrom0_xranges.append(chrom_xrange)
                    chrom0_source.append(source)

                elif genome_table.iloc[x, 1] == 1:

                    chrom1_xranges.append(chrom_xrange)
                    chrom1_source.append(source)    

            chrom0_color_list = []
            for x in chrom0_source:

                if x == 1:

                    chrom0_color_list.append("r")

                if x == 2:

                    chrom0_color_list.append("b")

            chrom1_color_list = []
            for x in chrom1_source:

                if x == 1:

                    chrom1_color_list.append("r")

                if x == 2:

                    chrom1_color_list.append("b")


            ax.broken_barh(chrom0_xranges, 
                          yrange = (ind_genome, ind_genome + 0.5),
                          facecolors = chrom0_color_list)

            ax.broken_barh(chrom1_xranges, 
                          yrange = (ind_genome + 0.5, ind_genome + 1),
                          facecolors = chrom1_color_list)
            
            
    def vis_single_genome(self, genome_index):

        figure, ax = plt.subplots(figsize = (10, 1))
        ax.set_xlim(0, self.seq_length)
        ax.set_ylim(0, 1)
        plt.grid(visible = False)
        
        plt.tick_params(left = False, right = False , labelleft = False)

        genome_table = self.genome_df_list[genome_index]

        mintick = np.arange(0, 1, 0.5)
        ax.set_yticks(mintick, minor=True)

        chrom0_xranges = []
        chrom1_xranges = []
        chrom0_source = []
        chrom1_source = []

        for x in range(len(genome_table)):

            chrom_xrange = (genome_table.iloc[x, 2], 
                            genome_table.iloc[x, 3] - genome_table.iloc[x, 2])
            source = genome_table.iloc[x, 4]

            if genome_table.iloc[x, 1] == 0:

                chrom0_xranges.append(chrom_xrange)
                chrom0_source.append(source)

            elif genome_table.iloc[x, 1] == 1:

                chrom1_xranges.append(chrom_xrange)
                chrom1_source.append(source)    
                
        chrom0_color_list = []
        for x in chrom0_source:

            if x == 1:

                chrom0_color_list.append("r")

            if x == 2:

                chrom0_color_list.append("b")

        chrom1_color_list = []
        for x in chrom1_source:

            if x == 1:

                chrom1_color_list.append("r")

            if x == 2:

                chrom1_color_list.append("b")


        ax.broken_barh(chrom0_xranges, 
                        yrange = (0, 0.5),
                        facecolors = chrom0_color_list)

        ax.broken_barh(chrom1_xranges, 
                        yrange = (0.5, 1),
                        facecolors = chrom1_color_list)
            
            
    def write_genome_df_list(self, 
                          output_filename):

        genome_df = pd.concat(self.genome_df_list)
        genome_df.to_csv(output_filename)

                          
sim1 = Sim(input_filename = "written_parent_frame.csv",
            seq_length = 100,
            recomb_rate = 0.01
         )

sim1.execute_sim()
#sim1.show_tree()
sim1.extract_genome_df_list()
sim1.vis_single_genome(0)
sim1.vis_genome_df_list()
sim1.write_sim(output_filename = "test ts output")
