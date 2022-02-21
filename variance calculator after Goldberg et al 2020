#to test variance directly
import numpy as np
import matplotlib.pyplot as plt


class VarCalc:
    
    def __init__(self, g, 
                s1f_0, s1m_0, s2f_0, s2m_0, 
                s1f, s1m, s2f, s2m,
                c11_0, c22_0, 
                c11, c1h, c12,
                ch1, chh, ch2,
                c21, c2h, c22):

        c12_0 = -c11_0
        c21_0 = -c22_0

        hf = 1 - s1f - s2f
        hm = 1 - s1m - s2m

        s1_0 = (s1f_0 + s1m_0) / 2
        s1 = (s1f + s1m) / 2
        s2 = (s2f + s2m) / 2
        h = (hf + hm) / 2
        
        #c22 doesn't appear in functions?!

        if c11 + c1h + c12 != 0 or\
        ch1 + chh + ch2 != 0 or\
        c21 + c2h + c22 != 0 or\
        c11 + ch1 + c21 != 0 or\
        c1h + chh + c2h != 0 or\
        c21 + c2h + c22 != 0:
            raise Exception('c values do not sum to zero!')

        predicted_var = []
        predicted_var.append(None)
        E_H_list = []
        E_H_list.append(None)

        for g in range(1, g + 1):
            if g == 1:
                #from eq 16
                E_H_list.append((s1f_0 * s1m_0) + c11_0 + (((s1f_0 * s2m_0) + (s2f_0 * s1m_0) + c12_0 + c21_0) / 2))
                
                #from eq 20
                predicted_var.append(((s1f_0 * (1 - s1f_0)) + (s1m_0 * (1 - s1m_0)) + (2 * c11_0)) / 4)

            elif g > 1:
                #from eq 17
                E_H_list.append(((s1f * s1m) + c11) +
                               (((s1f * s2m) + (s2f * s1m) + c12 + c21 + c1h + ch1) / 2) +
                               ((((s1f * hm) + (hf * s1m) + c1h + ch1) / 2) * (1 + E_H_list[g-1])) +
                               ((((hf * hm) + chh) / 2) * (E_H_list[g-1] + E_H_list[g-1])) +
                               ((((s2f * hm) + (hf * s2m) + ch2 + c2h) / 2) * E_H_list[g-1])
                               )                     
                
                #from eq 21
                predicted_var.append((((s1f * (1 - s1f)) + (s1m * (1 - s1m)) + (2 * c11)) / 4) +
                                    (((c1h + ch1 - (s1f * hf) - (s1m * hm)) / 2) * E_H_list[g-1]) +
                                    ((((hf * (1 - hf)) + (hm * (1 - hm)) + (2 * chh)) / 4) * (E_H_list[g-1] ** 2)) +
                                    (((hf + hm) / 4) * predicted_var[g-1])
                                    )

        x = range(len(predicted_var))
        var_predx.plot(predicted_var, color = ('red'))
        var_predx.set_ylim(0, .25)
        var_predx.set_xlim(0, len(x))
        plt.xticks(np.arange(min(x), max(x) + 1, 5))
        print(predicted_var[-1])

                  
predicted_varfig = plt.figure()
var_predx = predicted_varfig.add_subplot(111)
var_predx.set(xlabel = 'generations', ylabel = 'variance of S1 ancestry', title = 'Variance in Ancestry')

varcalc1 = VarCalc(g=20, 
                s1f_0=0.5, s1m_0=0.5, s2f_0=0.5, s2m_0=0.5, 
                s1f=0.2, s1m=0.2, s2f=0.2, s2m=0.2,
                c11_0=0, c22_0=0, 
                c11=0.02, c1h=-0.02, c12=0,
                ch1=-0.02, chh=0.02, ch2=0,
                c21=0, c2h=0, c22=0)

varcalc1 = VarCalc(g=20, 
                s1f_0=0.5, s1m_0=0.5, s2f_0=0.5, s2m_0=0.5, 
                s1f=0.2, s1m=0.2, s2f=0.2, s2m=0.2,
                c11_0=0, c22_0=0, 
                c11=0.02, c1h=-0.04, c12=0.02,
                ch1=-0.04, chh=0.02, ch2=0.02,
                c21=0.02, c2h=0.02, c22=-0.04)

varcalc1 = VarCalc(g=20, 
                s1f_0=0.5, s1m_0=0.5, s2f_0=0.5, s2m_0=0.5, 
                s1f=0.2, s1m=0.2, s2f=0.2, s2m=0.2,
                c11_0=0, c22_0=0, 
                c11=0.02, c1h=0.01, c12=-0.03,
                ch1=0.01, chh=0.02, ch2=-0.03,
                c21=-0.03, c2h=-0.03, c22=0.06)
