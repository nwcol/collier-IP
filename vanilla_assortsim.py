import pandas
import vanilla_coalescence_sim
import inference

seeds = []
gs = []
s1s = []
s2s = []
for x in range(10):

    sim = vanilla_coalescence_sim.VanillaCoalescenceSim(
        admix_time = 20, 
        pop_size = 200, 
        sample_size = 20, 
        cM_length = 300, 
        seq_length = 300000000, 
        s1_init = 0.5, 
        s2_init = 0.5)
    sim.execute_sim(output_dir = "vanillachrom_output/")
    print("simulation complete")

    inf = inference.TwoParamInference(input_dir = "vanillachrom_output/")
    mig, g, s1, s2 = inf.make_inference()
    gs.append(g)
    s1s.append(s1)
    s2s.append(s2)
    
dic = {"g": gs} #,
       #"s1": s1s,
       #"s2": s2s}
df = pandas.DataFrame(dic)
print(df)