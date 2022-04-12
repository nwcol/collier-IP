import pandas

import pop_sim
import coalescence_sim
import inference

trials = 2
sample_size = 40

seeds = []
gs = []
s1s = []
s2s = []

for x in range(trials):
    model = pop_sim.PopSim(
            g = 10, 
            N_init = 2000,
            N_const = 0,
            c11_0 = 0)
    seed = model.seed
    model.execute_model()
    model.write_sample_df(sample_size, "pop_sim_output/output.csv")
    print("pop_model complete")
    
    sim_1 = coalescence_sim.CoalescenceSim(
        seed = seed,
        input_filename = "pop_sim_output/output.csv",
        seq_length = 300000000,
        cM_length = 300)
    sim_1.make_sim()
    print("coalescence_sim complete")
    sim_1.write_genomes(output_dir = "coalescence_sim_output/")
    
    inf = inference.TwoParamInference(input_dir = "coalescence_sim_output/")
    mig, g, s1, s2 = inf.make_inference()
    
    seeds.append(seed)
    gs.append(g)
    s1s.append(s1)
    s2s.append(s2)

dic = {"seeds": seeds,
       "g": gs} #,
       #"s1": s1s,
       #"s2": s2s}
df = pandas.DataFrame(dic)
print(df)