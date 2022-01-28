I wrote most of this before actually writing the simulation, so it may be unclear or have actual discrepencies with the script.


Individuals are represented as lists (vectors?) with a sex dimension (0 = female, 1 = male) and an S1 ancestry dimension, which ranges from 0 to 1 and represents the fraction of the individual's ancestors originating from S1.

List of Variables Used: gen_0 = a list representing the founding population (representing a mixture of S1 and S2 individuals)

gen_gminus1 = the population at generation g - 1

gen_g = the population at generation g

1f0, 1m0, etc = the numbers of individuals from sources S1 and S2 that contribute to generation 0. There are no variables for h yet because it does not yet exist.

c11, c1h, c12 etc = cij is the bias in probability that a female of type i (1, h or 2) mates with a male of type j. All the cs in a row must sum to 0.

S1f, Hf, S2f etc = lists containing the individuals of a specific type and sex in g - 1. Used in selecting mates

s1f, hf, s2f etc = the fraction of indvidiuals of a specific type and sex in g - 1. Used to calculate Ps

P11, P1h, P12 etc = the probability Pij that a random individual from generation g had parents i (female) and j (male) from generation g - 1. Calculated from the fractions

selection = the mating which is selected in the mate function

offspring = a new member of generation g, which is constructed appropriately and then appended to gen_g

hybrid_father, hybrid_mother = used to select a random member of H to calculate the ancestry of an individual resulting from mating with or within H

g = the number of generations to go to

gens = a list containing all generations - index = generation eg. gens[0] is gen_0

ancestries = a list of all ancestry fractions from a generation

Method: to produce generation g, a number of random mating events equal to the n individuals in generation g - 1 (to produce an equally-sized population--- each event produces one offspring in generaiton g) are processed. Essentially, for each x in range(n), a weighted choice is made based upon probabilities based on the composition of g - 1. Then, two individuals from appropriate populations are selected. They produce one offspring of random sex, whose ancestry number is the average of the parents (eg. parents of ancestry 1 and .5 produce a child of ancestry .75). Finally this offspring is appended into generation g.
