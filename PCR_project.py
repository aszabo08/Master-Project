"""
This application visualizes primer binding and misbinding to a substrate, creating substrate-primer products.
"S" stands for the concentration of the substrate, "P" and "SP" for the concentration of the primer and the bound state of the substrate and primer.
The concentrations for all the above mentioned species are given in mole/liter.
The k parameters (k1, k2, k3) are rate constants in binding and misbinding processes.
"""


from scipy.integrate import odeint

import matplotlib.pyplot as plt

import numpy as np



species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]


values = [0 for i in range(17)]


values[0] = 3.0769e-2       # concentration of plasmid (S1S2) in uM
values[1] = 0               # concentration of S1 in uM
values[2] = 0               # concentration of S2 in uM

values[3] = 2             # concentration of P1 in uM
values[4] = 2            # concentration of P2 in uM

values[5] = 0               # concentration of S1P2
values[6] = 0               # concentration of S2P1
values[7] = 0.2             # concentration of E in uM
values[8] = 0               # concentration of S1P2E
values[9] = 0               # concentration of S2P1E
values[10] = 200           # concentration of dNTP in uM
values[11] = 0              # concentration of Q1
values[12] = 0              # concentration of Q2
values[13] = 0              # concentration of S1Q2E
values[14] = 0              # concentration of S2Q1E
values[15] = 0              # concentration of S2Q1
values[16] = 0              # concentration of S1Q2


initial_dNTP = values[10]

Tden = 369.15               # 96 degree

Tanneal = 303.15            # 30 degree


Text = 343.15               # 70 degree

tden = 10                   # seconds
tanneal = 10                # seconds
text = 10                   # seconds


total = tden + tanneal + text

number_cycles = 1

steps = 1


amplicon_length = 1000
primer_length = 15

n = 10

extended_primer = primer_length + n

extended_length = amplicon_length - extended_primer

time = np.linspace(0, number_cycles * (tden + tanneal + text), number_cycles * (tden + tanneal + text) * steps)  # for every second "steps" points are distinguished

number_time_points = total * steps * number_cycles



# 1 entropy unit = 4.184 J/ K m


Tm_S1S2 = 358.15                # 85 degrees

Tm_primer = 308.15              # 35 degrees

Tm_extended_primer = 313.15     # 40 degree


#dS = -1256.0395 # j/m

dS = - 52.7184e-2





R = 8.314e-3        # Gas contant in  J / K mol

#kB = 1.38064852e-23  # Boltzmann constant

def number_of_molecules(x):

    """This function calculates the number of molecules with given concentration in uM in a 50 ul tube"""

    umol_per_ul = x * 1e-6     # 1 uM = 10e-6 mol/l    which is equal to 10e-6 umol/ul

    umol = umol_per_ul * 50       # number of umol in a 50 ul tube

    molecules = umol * 6.022140857 * 1e+17     # in 1 mol there is 6.022140857 * 10e+23 molecules  -> in 1 umol there is 6.022140857 * 10e+23 * 10e-6 = 6.022140857 * 10e+17

    return molecules


def total_molecule_length(x, y):

    """This function calculates the summary of the lengths of every molecules within a species in nucleotides"""

    return number_of_molecules(x) * y


bp_S1S2 = total_molecule_length(values[0], amplicon_length)


#all_nucleotide = 2 * total_molecule_length(values[0], amplicon_length) + total_molecule_length(values[1], amplicon_length) + total_molecule_length(values[2], amplicon_length) + total_molecule_length(values[3], primer_length) total_molecule_length(values[4], primer_length)          # in the tube


def all_nucleotide():

    length = [2 * amplicon_length, amplicon_length, amplicon_length, primer_length, primer_length, amplicon_length + primer_length, amplicon_length + primer_length, 0, amplicon_length + primer_length, amplicon_length + primer_length, 1, extended_primer, extended_primer, amplicon_length + extended_primer,  amplicon_length + extended_primer, amplicon_length + extended_primer, amplicon_length + extended_primer]

    total_number = 0

    for i in range(len(values)):


        total_number = total_number + total_molecule_length(values[i], length[i])



    return total_number


def used_dNTP(x,y,z):


    dNTP_in_primer_ext1_s1 = number_of_molecules(x) * n

    dNTP_in_primer_ext1_s2 = number_of_molecules(y) * n

    dNTP_in_primer_ext2_s1s2 = number_of_molecules(z) * extended_length


    used_dNTP = dNTP_in_primer_ext1_s1 + dNTP_in_primer_ext1_s2 + dNTP_in_primer_ext2_s1s2


    return used_dNTP


def taq_nt_per_s(temperature):



    first_interval = np.linspace(0, 0.25, (22 - 0 + 1))
    second_interval = np.linspace(0.25, 1.5, (37 - 22 + 1))
    third_interval = np.linspace(1.5, 24, (55 - 37 + 1))
    forth_interval = np.linspace(24, 60, (70 - 55 + 1))
    fifth_interval = np.linspace(60, 150, (75 - 70 + 1))
    sixth_interval = np.full((1, (79 - 75)), 150)
    seventh_interval = np.linspace(150, 0, (90 - 80 + 1))



    interval_values = np.concatenate((first_interval, second_interval, third_interval, forth_interval, fifth_interval), axis=None)


    nt_per_s = np.concatenate((np.unique(interval_values), sixth_interval, seventh_interval) , axis=None)


    taq_temperature = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(90), (90 - 0 + 1))

    taq_index = list(taq_temperature).index(temperature)

    return nt_per_s[taq_index]



def denaturation(values, t, T, dGs):


    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]

    #print("S1S2, S1, S2", S1S2, S1, S2 )


    kf1 = 1


    exponent_1 = dGs[0]/(R*T)

    exponent_1 = np.clip(exponent_1, a_min= None,  a_max= 29.9336 )

    kr1 = kf1 * np.exp(exponent_1)

    #print("kr1 origin", kr1)

    kr1 = np.clip(kr1, a_min= -1e+14, a_max= 1e+14)

    if np.abs(kr1<1e-14):
        kr1 = 0

    #print("kr1 second", kr1)



    rate_den = - kr1 * S1S2 + kf1 * S1 * S2

    #print("rateden first", rate_den)

    rate_den = np.clip(rate_den,  a_min= -1e+14, a_max= 1e+14 )

    #print("rateden", rate_den)

    y = np.zeros(17)

    y[0] = rate_den

    y[1] = -rate_den

    y[2] = -rate_den



    return y



def primer_binding_1(values, t, T, dGs):


    """
    Describes the differential equations for S, P and SP in a dynamic biological reaction
    Args:
        values: value of S, P and SP after t time period
        t: time period
        k1: rate constant for binding ( S + P -> SP)
        k2: rate constant for dissociation ( SP -> S + P)
    Returns:
        An array with the result of the differential equations for S, P and SP over t time
    """


    S1 = values[1]
    S2 = values[2]
    P1 = values[3]
    P2 = values[4]
    S1P2 = values[5]
    S2P1 = values[6]


    kf2 = 1



    exponent_2a = dGs[1]/(R*T)

    #print("exponent_2a origin", exponent_2a)

    exponent_2a = np.clip(exponent_2a, a_min= None, a_max= 29.9336)

    #print("exponent_2a second", exponent_2a)

    kr2a = kf2 * np.exp(exponent_2a)

    kr2a = np.clip(kr2a,  a_min= -1e+14, a_max= 1e+14 )

    if np.abs(kr2a < 1e-14):
        kr2a = 0



    #print("kr2a second", kr2a)



    rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2


    rate_S1P2_bind = np.clip(rate_S1P2_bind,  a_min= -1e+14, a_max= 1e+14 )

    #print("rate S1P1 bind", rate_S1P2_bind)

    exponent_2b = dGs[1]/(R*T)

    #print("exponent_2b origin", exponent_2b)

    exponent_2b = np.clip(exponent_2b, a_min= None, a_max= 29.9336)

    #print("exponent_2b second", exponent_2b)

    kr2b = kf2 * np.exp(exponent_2b)


    kr2b = np.clip(kr2b,  a_min= -1e+14, a_max= 1e+14 )


    if np.abs(kr2b < 1e-14):
        kr2b = 0


    #print("kr2b", kr2b)

    rate_S2P1_bind = kf2 * S2 * P1 - kr2b * S2P1


    rate_S2P1_bind = np.clip(rate_S2P1_bind,  a_min= -1e+14, a_max= 1e+14 )

    #print("rate_s2P1", rate_S2P1_bind)


    y = np.zeros(17)


    y[1] = - rate_S1P2_bind
    y[2] = - rate_S2P1_bind
    y[3] = - rate_S2P1_bind
    y[4] = - rate_S1P2_bind
    y[5] = rate_S1P2_bind
    y[6] = rate_S2P1_bind


    return y



def primer_binding_2(values, t, T, dGs):


    """
    Describes the differential equations for S, P and SP in a dynamic biological reaction
    Args:
        values: value of S, P and SP after t time period
        t: time period
        k1: rate constant for binding ( S + P -> SP)
        k2: rate constant for dissociation ( SP -> S + P)
    Returns:
        An array with the result of the differential equations for S, P and SP over t time
    """


    S1 = values[1]
    S2 = values[2]
    Q1 = values[11]
    Q2 = values[12]
    S1Q2 = values[15]
    S2Q1 = values[16]


    kf5 = 1

    exponent_5a = dGs[2]/(R*T)

    exponent_5a = np.clip(exponent_5a, a_min= None, a_max= 29.9336)

    #print("exponent_5a", exponent_5a)

    kr5a = kf5 * np.exp(exponent_5a)

    kr5a = np.clip(kr5a,  a_min= -1e+14, a_max= 1e+14 )

    if np.abs(kr5a < 1e-14):
        kr5a = 0

    #print("kr5a", kr5a)

    rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    rate_S1Q2_bind = np.clip(rate_S1Q2_bind,  a_min= -1e+14, a_max= 1e+14 )

    #print("rate_S1Q2_bind", rate_S1Q2_bind)

    exponent_5b = dGs[2]/(R*T)

    exponent_5b = np.clip(exponent_5b, a_min= None, a_max= 29.9336)


    #print("exponent_5b", exponent_5b)

    kr5b = kf5 * np.exp(exponent_5b)


    kr5b = np.clip(kr5b,  a_min= -1e+14, a_max= 1e+14 )


    if np.abs(kr5b < 1e-14):
        kr5b = 0

    #print("kr5b", kr5b)

    rate_S2Q1_bind = kf5 * S2 * Q1 - kr5b * S2Q1

    rate_S2Q1_bind = np.clip(rate_S2Q1_bind,  a_min= -1e+14, a_max= 1e+14 )

    #print("rate_S2Q1_bind", rate_S2Q1_bind)

    y = np.zeros(17)


    y[1] = - rate_S1Q2_bind
    y[2] = - rate_S2Q1_bind
    y[11] = - rate_S2Q1_bind
    y[12] = - rate_S1Q2_bind
    y[15] = rate_S1Q2_bind
    y[16] = rate_S2Q1_bind


    return y




def polymerase_binding_1(values, t, T, dGs):


    S1P2 = values[5]
    S2P1 = values[6]
    E = values[7]
    S1P2E = values[8]
    S2P1E = values[9]


    kf3 = 1

    exponent_3a = dGs[3]/(R*T)

    exponent_3a = np.clip(exponent_3a, a_min= None, a_max= 29.9336)

    kr3a = kf3 * np.exp(exponent_3a)

    kr3a = np.clip(kr3a,  a_min= -1e+14, a_max= 1e+14 )

    if np.abs(kr3a<1e-14):
        kr3a = 0

    rate_poly_S1P2_bind = kf3 * S1P2 * E - kr3a * S1P2E

    rate_poly_S1P2_bind = np.clip(rate_poly_S1P2_bind,  a_min= -1e+14, a_max= 1e+14 )

    exponent_3b = dGs[3]/(R*T)

    exponent_3b = np.clip(exponent_3b, a_min= None, a_max= 29.9336)


    kr3b = kf3 * np.exp(exponent_3b)

    kr3b = np.clip(kr3b,  a_min= -1e+14, a_max= 1e+14 )

    if np.abs(kr3b<1e-14):
        kr3b = 0

    rate_poly_S2P1_bind = kf3 * S2P1 * E - kr3b * S2P1E

    rate_poly_S2P1_bind = np.clip(rate_poly_S2P1_bind,  a_min= -1e+14, a_max= 1e+14)

    enzyme_binding = - rate_poly_S1P2_bind - rate_poly_S2P1_bind

    enzyme_binding = np.clip(enzyme_binding,  a_min= -1e+14, a_max= 1e+14 )


    y = np.zeros(17)

    y[5] = - rate_poly_S1P2_bind
    y[6] = - rate_poly_S2P1_bind
    y[7] = enzyme_binding
    y[8] = rate_poly_S1P2_bind
    y[9] = rate_poly_S2P1_bind


    return y



def polymerase_binding_2(values, t, T, dGs):


    S1Q2 = values[15]
    S2Q1 = values[16]
    E = values[7]
    S1Q2E = values[13]
    S2Q1E = values[14]


    kf4 = 1

    exponent_4a = dGs[3]/(R*T)

    exponent_4a = np.clip(exponent_4a, a_min= None, a_max= 29.9336)


    kr4a = kf4 * np.exp(exponent_4a)

    kr4a = np.clip(kr4a,  a_min= -1e+14, a_max= 1e+14)

    if np.abs(kr4a<1e-14):
        kr4a = 0

    rate_poly_S1Q2_bind = kf4 * S1Q2 * E - kr4a * S1Q2E

    rate_poly_S1Q2_bind = np.clip(rate_poly_S1Q2_bind,  a_min= -1e+14, a_max= 1e+14)



    exponent_4b = dGs[3]/(R*T)

    exponent_4b = np.clip(exponent_4b, a_min= None, a_max= 29.9336)

    #print("exponent_4b", exponent_4b)


    kr4b = kf4 * np.exp(exponent_4b)

    kr4b = np.clip(kr4b,  a_min= -1e+14, a_max= 1e+14 )

    if np.abs(kr4b<1e-14):
        kr4b = 0

    rate_poly_S2Q1_bind = kf4 * S2Q1 * E - kr4b * S2Q1E

    rate_poly_S2Q1_bind = np.clip(rate_poly_S2Q1_bind,  a_min= -1e+14, a_max= 1e+14)

    enzyme = - rate_poly_S1Q2_bind - rate_poly_S2Q1_bind

    enzyme = np.clip(enzyme, a_min= -1e+14, a_max= 1e+14 )


    y = np.zeros(17)

    y[15] = - rate_poly_S1Q2_bind
    y[16] = - rate_poly_S2Q1_bind
    y[7] = enzyme
    y[13] = rate_poly_S1Q2_bind
    y[14] = rate_poly_S2Q1_bind


    return y







def primer_ext_1(values, t, T, dGs):

    """ The primer is extended by a few nucleotides to ensure the primer binding
        to the substrate without dissociation when reaching the extension temperature

        In this model 99.5 % of the primer - substrate complexes stay together, while 0.5 % of them will melt
        at the extension temperature"""


    S1P2E = values[8]
    S2P1E = values[9]
    dNTP = values[10]

    #   S1P2E + n * dNTP -> S1Q2E
    #   S2P1E + n * dNTP -> S2Q1E

    ce = 1000   # [dNTP/s] concentration of polymerase enzyme

    rate_ext_1 = (ce / n) * S1P2E * dNTP

    rate_ext_1 = np.clip(rate_ext_1,  a_min= -1e+14, a_max= 1e+14 )

    #print("rate_exp_1", rate_ext_1)

    rate_ext_2 = ce / n * S2P1E * dNTP

    rate_ext_2 = np.clip(rate_ext_2,  a_min= -1e+14, a_max= 1e+14)


    #print("rate_ext_2", rate_ext_2)

    nucleotide = - n * rate_ext_1 - n * rate_ext_2

    nucleotide = np.clip(nucleotide,  a_min= -1e+14, a_max= 1e+14)


    primer_ext_1.counter += 1


    y = np.zeros(17)


    y[8] = - rate_ext_1    # concentration of S1P2E
    y[9] = - rate_ext_2   # concentration of S2P1E
    y[10] = nucleotide   # concentration of dNTP
    y[13] = rate_ext_1    # concentration of S1Q2E
    y[14] = rate_ext_2   # concentration of S2Q1E


    return y


def primer_ext_2(values, t, T, dGs):



    S1Q2E = values[13]
    S2Q1E = values[14]
    dNTP = values[10]

    ce_Q = 1000         # enzyme concentration

    # reaction: S1Q2E + extended_length * dNTP ---> S1S2 + E

    rate_ext_Q1 = (ce_Q / extended_length) * S1Q2E * dNTP

    rate_ext_Q1 = np.clip(rate_ext_Q1,  a_min= -1e+14, a_max= 1e+14)


    # reaction: S2Q1E + extended_length * dNTP ---> S1S2 + E

    rate_ext_Q2 = (ce_Q / extended_length) * S2Q1E * dNTP

    rate_ext_Q2 = np.clip(rate_ext_Q2,  a_min= -1e+14, a_max= 1e+14)


    nucleotide_Q = - extended_length * rate_ext_Q1 - extended_length * rate_ext_Q2

    nucleotide_Q = np.clip(nucleotide_Q,  a_min= -1e+14, a_max= 1e+14)

    primer_ext_2.counter += 1

    y = np.zeros(17)

    y[0] = rate_ext_Q1 + rate_ext_Q2
    y[7] = rate_ext_Q1 + rate_ext_Q2
    y[10] = nucleotide_Q                            # concentration of dNTP
    y[13] = - rate_ext_Q1                           # concentration of S1Q2E
    y[14] = - rate_ext_Q2                           # concentration of S2Q1E




    return y


def PCR_reaction(values, t, T, dGs):

    summary = denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + primer_binding_2(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + primer_ext_2(values, t, T, dGs)

    summary = np.clip(summary,  a_min= -1e+14, a_max= 1e+14)

    summary[np.abs(summary) < 1e-14] = 0

    return summary

#
def iterative_PCR(values, t, T, dGs):


    summary = 0

    #result = np.zeros((7, 17))


    for i in range(len(functions)):

        summary = summary + functions[i](values, t, T, dGs)



    return summary



def new_dntp():

    used_dntp = primer_ext_1.counter * 2 * n + primer_ext_2.counter * 2 * extended_length

    return used_dntp



if __name__ == '__main__':

    functions = [denaturation, primer_binding_1, polymerase_binding_1, primer_ext_1, primer_ext_2, polymerase_binding_2, primer_binding_2]

    print("all nuc", all_nucleotide())

    first = all_nucleotide()

    concentration = np.empty((number_time_points, 17))

    dGs = [0,0,0,0]

    primer_ext_1.counter = 0

    primer_ext_2.counter = 0

    dNTP_check = 0

    print(number_of_molecules(initial_dNTP))

    result = np.empty(0)



    for i in range(number_cycles):


        dGs[0] = (Tm_S1S2 - Tden) * dS

        dGs[1] = (Tm_primer - Tden) * dS

        dGs[2] = (Tm_extended_primer - Tden) * dS

        dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tden * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+14)

        print("dGs_den", dGs)

        before_level = [values[13], values[14], values[0]]



        integration_den = odeint(PCR_reaction, values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs))

        concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den

        dGs[0] = (Tm_S1S2 - Tanneal) * dS

        dGs[1] = (Tm_primer - Tanneal) * dS

        dGs[2] = (Tm_extended_primer - Tanneal) * dS

        dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tanneal * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+14)

        print("dGs_anneals", dGs)

        integration_anneal = odeint(PCR_reaction, integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs))

        concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal

        dGs[0] = (Tm_S1S2 - Text) * dS

        dGs[1] = (Tm_primer - Text) * dS

        dGs[2] = (Tm_extended_primer - Text) * dS

        dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Text * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+14)

        print("dGs_text", dGs)

        integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs))

        concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


        values = integration_ext[-1]

        after_level = [values[13], values[14], values[0]]

        #dNTP_check = dNTP_check + used_dNTP(after_level[0] - before_level[0],after_level[1] - before_level[1], after_level[2] - before_level[2])

        #print("dntp check", number_of_molecules(initial_dNTP) - dNTP_check)

        print("loop new dntp", new_dntp())

        dNTP_check = dNTP_check + new_dntp()

        print("dntp check", dNTP_check)

        final_level = number_of_molecules(initial_dNTP) - dNTP_check

        print("final level", final_level)

        print("dntp value", number_of_molecules(values[10]))



        print(all_nucleotide())


    print("newdntp", new_dntp())



    print("The concentration of the 17 species at the end of the process:", values)

    print("counter", primer_ext_1.counter)

    print("counter2", primer_ext_2.counter)

    print(total_molecule_length(1,10))

    print("all nuc", all_nucleotide())

    second = all_nucleotide()

    print("diff", first - second)


    #print("concentration", concentration)


    #Plotting the concentrartions over time

    plt.figure(1)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)

    plots = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]


    for i in range(10):

        plt.subplot(2, 5, i+1)


        plt.plot(time, concentration[:, plots[i]])


        plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time")
        plt.ylabel("Concentration")



    plt.show()















