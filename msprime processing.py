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

        pedigree_df = pandas.read_csv(self.input_filename)
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
                y_axis = True,
                y_ticks = [0, 5, 10, 15, 20, 25, 30]
            )
            show(SVG(sim_svg))
            
            
    def extract_genome_list(self):
        
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

        genome_list = []

        for ind in range(len(sample_nodes)):
            
            ind_nodes = sample_nodes[ind]

            chroms = []
            starts = []
            ends = []
            sources = []

            #search for each chromosome     
            for chrom in range(2):

                for row in range(len(self.sim_ts.tables.edges)):

                    selected_row = self.sim_ts.tables.edges[row]

                    if selected_row.child == ind_nodes[chrom]:

                        chroms.append(chrom)
                        starts.append(selected_row.left)
                        ends.append(selected_row.right)
                        parent_node = selected_row.parent

                        source = None
                        #temporary- to prevent infinite recursion if something breaks
                        counter = 0

                        while source == None and counter <20:

                            counter = counter + 1

                            child_node = parent_node

                            for y in range(len(self.sim_ts.tables.edges)):

                                selected_row = self.sim_ts.tables.edges[y]

                                if selected_row.child == child_node:

                                    parent_node = selected_row.parent

                            selected_row = self.sim_ts.tables.nodes[parent_node]
                            selected_ind = selected_row.individual
                            selected_row = self.sim_ts.tables.individuals[selected_ind]
                            selected_parents = selected_row.parents

                            if selected_parents[0] == -1:

                                source_id = self.sim_ts.tables.nodes[parent_node].individual
                                source_metadata = self.sim_ts.tables.individuals[source_id].metadata
                                source = int(source_metadata['$comment'])
                                sources.append(source)

                            else:

                                pass

                        source = None
                        
            ind_id = [sample_nodes[ind]] * len(chroms)

            ind_genome = {"chromosome": chroms,
                          "start": starts,
                          "end": ends,
                          "source": sources,
                          "id": ind_id
                        }

            ind_genome_df = pd.DataFrame(ind_genome)
            ind_genome_df = ind_genome_df.sort_values(by = ["chromosome", "start"])
            genome_list.append(ind_genome_df)
        
        self.genome_list = genome_list

    
    def vis_genome_list(self):

        figure, ax = plt.subplots(figsize = (10, 0.25 * len(self.genome_list)))
        ax.set_ylim(0, len(self.genome_list))
        plt.grid(visible = True, axis = 'y', color = 'black', )

        for ind_genome in range(len(self.genome_list)):

            genome_table = self.genome_list[ind_genome]

            majtick = np.arange(0, len(self.genome_list), 1)
            mintick = np.arange(0, len(self.genome_list), 0.5)
            ax.set_yticks(majtick)
            ax.set_yticks(mintick, minor=True)

            labels = ['0', '1']

            chrom0_xranges = []
            chrom1_xranges = []
            chrom0_source = []
            chrom1_source = []

            for x in range(len(genome_table)):

                chrom_xrange = (genome_table.iloc[x, 1], 
                                genome_table.iloc[x, 2] - genome_table.iloc[x, 1])
                source = genome_table.iloc[x, 3]

                if genome_table.iloc[x, 0] == 0:

                    chrom0_xranges.append(chrom_xrange)
                    chrom0_source.append(source)

                elif genome_table.iloc[x, 0] == 1:

                    chrom1_xranges.append(chrom_xrange)
                    chrom1_source.append(source)    

            chrom0_color_list = []
            for x in chrom0_source:

                if x == 1:

                    chrom0_color_list.append("tab:red")

                if x == 2:

                    chrom0_color_list.append("tab:blue")

            chrom1_color_list = []
            for x in chrom1_source:

                if x == 1:

                    chrom1_color_list.append("tab:red")

                if x == 2:

                    chrom1_color_list.append("tab:blue")


            ax.broken_barh(chrom0_xranges, 
                          yrange = (ind_genome, ind_genome + 0.5),
                          facecolors = chrom0_color_list)

            ax.broken_barh(chrom1_xranges, 
                          yrange = (ind_genome + 0.5, ind_genome + 1),
                          facecolors = chrom1_color_list)
            
            
    def write_genome_list(self, 
                          output_filename):

        genome_df = pd.concat(self.genome_list)
        genome_df.to_csv(output_filename)
        return(genome_df)
                          
sim1 = Sim(input_filename = "written_parent_frame.csv",
            seq_length = 10,
            recomb_rate = 0.01
         )

sim1.execute_sim()
sim1.extract_genome_list()
sim1.vis_genome_list()   
sim1.write_genome_list(output_filename = "test of genome output.csv")
