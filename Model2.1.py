import random
import pandas as pd
import copy

def new_population(S1f_init, S1m_init, S2f_init, S2m_init):
    '''generates an initial generation 0 based on initial parameters'''
    
    gen_0 = []
    for x in range(S1f_init):
        gen_0.append([0, 1, 1])
    for x in range(S1m_init):
        gen_0.append([1, 1, 1])
    for x in range(S2f_init):
        gen_0.append([0, 0, 2])
    for x in range(S2m_init):
        gen_0.append([1, 0, 2])
        
    return(gen_0)
        
def constant_contribution(gen_g_minus1, S1f_const, S1m_const, S2f_const, S2m_const):
    '''appends constant contributions to gen_g_minus1'''
    
    for x in range(S1f_const):
        gen_g_minus1.append([0, 1, 1])
    for x in range(S1m_const):
        gen_g_minus1.append([1, 1, 1])
    for x in range(S2f_const):
        gen_g_minus1.append([0, 0, 2])
    for x in range(S2m_const):
        gen_g_minus1.append([1, 0, 2])
        
    return(gen_g_minus1)

def find_probs_gen1(gen_0):
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
    
    P11 = s1f * s1m + c11_0
    P12 = s1f * s2m + c12_0
    
    P21 = s2f * s1m + c21_0
    P22 = s2f * s2m + c22_0
    
    P1h, Ph1, Phh, Ph2, P2h = 0, 0, 0, 0, 0
    
    prob_list = [P11, P1h, P12, Ph1, Phh, Ph2, P21, P2h, P22]
    return(prob_list)

def find_probs(gen_g_minus1):
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
    hm = Hm/Nm, 
    s2m = S2m/Nm
    
    P11 = s1f * s1m + c11
    P1h = s1f * hm + c1h
    P12 = s1f * s2m + c12
    
    Ph1 = hf * s1m + ch1
    Phh = hf * hm + chh
    Ph2 = hf * s2m + ch2
    
    P21 = s2f * s1m + c21
    P2h = s2f * hm + c2h
    P22 = s2f * s2m + c22
    
    #interim way to remove negative probabilities
    prob_list = [P11, P1h, P12, Ph1, Phh, Ph2, P21, P2h, P22]
    #for x in range(9):
     #   if prob_list[x] < 0:
     #       prob_list[x] = 0          
    
    return(prob_list)

def mate(gen_g_minus1, is_gen_1):
    '''creates a new generation'''

    #define probabilities- gen 1 is a special case with its own prob function
    if is_gen_1 == True:
        prob_list = find_probs_gen1(gen_g_minus1)
    else:
        prob_list = find_probs(gen_g_minus1)
    
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
    
    for x in range(N):
        
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

def tabulate_parents_gen_0(gen_0):
    '''tabulates parentage for gen g individuals. x is generation number'''
    
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

def tabulate_parents(gen_g, gen_g_parents, x):
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

def generations(g):
    '''produces a list of generations'''
    gens = []
    parent_frame_list = []
    
    gen_0 = new_population(S1f_init, S1m_init, S2f_init, S2m_init)
    gens.append(gen_0) 
    gen_g_minus1 = copy.deepcopy(gen_0)
    gen_0_parent_frame = tabulate_parents_gen_0(gen_0)
    parent_frame_list.append(gen_0_parent_frame)
    
    for x in range(1, g + 1):
        
        if x == 1:
            gen_g, gen_g_parents = mate(gen_g_minus1, True)
            gens.append(gen_g)
            
        else:
            del(gen_g_minus1[(N - N_const):])
            gen_g_minus1 = constant_contribution(gen_g_minus1, S1f_const, S1m_const, S2f_const, S2m_const)
            gen_g, gen_g_parents = mate(gen_g_minus1, False)
            gens.append(gen_g)
            
        gen_g_parent_frame = tabulate_parents(gen_g, gen_g_parents, x)
        parent_frame_list.append(gen_g_parent_frame)
        
        gen_g_minus1 = copy.deepcopy(gen_g)
    
    parent_frames = pd.concat(parent_frame_list)
    
    return(gens, parent_frames)

def manygenerations(runs):
    
    manygens = []
    many_parent_frames_list = []
    
    for x in range(runs):
        gens, parent_frames = generations(g)
        manygens.append(gens)
        many_parent_frames_list.append(parent_frames)
        
    many_parent_frames = pd.concat(many_parent_frames_list) 
    
    return(manygens, many_parent_frames)
    
#Definitions of Constants:

#values for generation 0
c11_0, c12_0 = 0, 0
c21_0, c22_0 = 0, 0

#values for subsequent generations
c11, c1h, c12 = -0.04, 0.02, 0.02
ch1, chh, ch2 = -0.01, 0.02, -0.01
c21, c2h, c22 = -0.01, -0.01, 0.02

g = 20
runs = 100

S1f_init = 500
S1m_init = 500
S2f_init = 500
S2m_init = 500

S1f_const = 200
S1m_const = 200
S2f_const = 200
S2m_const = 200

N_const = S1f_const + S1m_const + S2f_const + S2m_const
N = S1f_init + S1m_init + S2f_init + S2m_init
Nf = S1f_init + S2f_init
Nm = S1m_init + S2m_init
N_const = S1f_const + S1m_const + S2f_const + S2m_const

S1 = (S1f_const + S1m_const) / N
S2 = (S2f_const + S2f_const) / N

manygens, many_parent_frames = manygenerations(runs)

%store S1
%store S2

%store g
%store manygens
%store N
%store N_const
%store c11_0
%store c11
%store many_parent_frames
