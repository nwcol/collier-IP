import numpy as np
import pandas
import tskit
import msprime
from IPython.display import SVG

pedigree_df = pandas.read_csv('written_parent_frame.csv')
pandas.set_option("display.max_rows", None, "display.max_columns", None)

generic_female_1 = -1
generic_male_1 = -1
generic_female_2 = -1
generic_male_2 = -1
max_gen = pedigree_df.iloc[len(pedigree_df) - 1, 1]

my_pedigree = msprime.PedigreeBuilder()

for x in range(len(pedigree_df)):
    female_parent_id = pedigree_df.iloc[x, 5]
    male_parent_id = pedigree_df.iloc[x, 6]
    
    if female_parent_id == -1:
        female_parent = generic_female_1
    elif female_parent_id == -2:
        female_parent = generic_female_2
    else:
        female_parent = female_parent_id
    
    if male_parent_id == -1:
        male_parent = generic_male_1
    elif male_parent_id == -2:
        male_parent = generic_male_2
    else:
        male_parent = male_parent_id
        
    inst_parents = [female_parent, male_parent]
    inst_time = max_gen - pedigree_df.iloc[x, 1]
    my_pedigree.add_individual(time=inst_time, parents=inst_parents, population=0)
    
my_pedigree.tables.individuals
