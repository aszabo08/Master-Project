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


#values[0] = 3.0769e-2       # concentration of plasmid (S1S2) in uM
values[0] = 3.0769e-2
values[1] = 0               # concentration of S1 in uM
values[2] = 0               # concentration of S2 in uM



#values[3] = 2             # concentration of P1 in uM
#values[4] = 2            # concentration of P2 in uM

values[3] = 8            # concentration of P1 in uM
values[4] = 8

values[5] = 0               # concentration of S1P2
values[6] = 0               # concentration of S2P1
values[7] = 0.2             # concentration of E in uM
values[8] = 0               # concentration of S1P2E
values[9] = 0               # concentration of S2P1E
#values[10] = 200           # concentration of dNTP in uM

values[10] = 10000



values[11] = 0              # concentration of Q1
values[12] = 0              # concentration of Q2
values[13] = 0              # concentration of S1Q2E
values[14] = 0              # concentration of S2Q1E
values[15] = 0              # concentration of S2Q1
values[16] = 0              # concentration of S1Q2


initial_dNTP = values[10]

initial_fixed = values

functions_name = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer_ext_2"]

#Tden = 343.15               # 96 degree

Tden = 369.15

Tanneal = 303.15           # 30 degree
#Tanneal = 323.15

Text = 343.15               # 70 degree

#Text = 323.15


tden = 10                   # seconds
tanneal = 10                # seconds
text = 10                   # seconds


total = tden + tanneal + text

number_cycles = 30

steps = 1


amplicon_length = 1000
primer_length = 15

n = 10

extended_primer = primer_length + n

extended_length = amplicon_length - extended_primer

time = np.linspace(0, number_cycles * (tden + tanneal + text), number_cycles * (tden + tanneal + text) * steps)  # for every second "steps" points are distinguished

number_time_points = total * steps * number_cycles


# 1 entropy unit = 4.184 J/ K m


Tm_S1S2 = 363.15                # 90 degrees

Tm_primer = 308.15              # 35 degrees

Tm_extended_primer = 313.15     # 40 degree

Tm_enzyme = 353.15              # 80 degree


#dS = -1256.0395 # j/m

#dS = - 52.7184e-2


#dS = - 52.7184e-2/2         # -0.263592

#dS = - 52.7184e-1


dS = - 2


max_exponent = 20       # e on the power of 15 = 3,269,017.3724702

min_clip = -1e+30

max_clip = 1e+30



forward_rate = 1e-5


                        






R = 8.314e-3        # Gas contant in  J / K mol

#kB = 1.38064852e-23  # Boltzmann constant



def celsius_to_Kelvin(x):

    return x + 273.15



def clipping(kf, kr):

    return (kf / kr) * np.clip(kr, a_min= min_clip, a_max= max_clip), np.clip(kr, a_min= min_clip, a_max= max_clip)


def rate_clipping(x):

    if x != np.clip(x, a_min= min_clip, a_max= max_clip):

        print("Clipped rate value:", x)


    return np.clip(x, a_min= min_clip, a_max= max_clip)


def exponent_clipping(x):

    if x != np.clip(x, a_min= None, a_max= max_exponent):

        print("Clipped exponent value:", x)


    return np.clip(x, a_min= None, a_max= max_exponent)









def number_of_molecules(x):

    """This function calculates the number of molecules with given concentration in uM in a 50 ul tube"""

    umol_per_ul = x * 1e-6     # 1 uM = 1e-6 mol/l    which is equal to 1e-6 umol/ul

    umol = umol_per_ul * 50       # number of umol in a 50 ul tube

    molecules = umol * 6.022140857 * 1e+17     # in 1 mol there is 6.022140857 * 1e+23 molecules  -> in 1 umol there is 6.022140857 * 1e+23 * 10e-6 = 6.022140857 * 1e+17

    return molecules


def total_molecule_length(x, y):

    """This function calculates the summary of the lengths of every molecules within a species in nucleotides"""

    return number_of_molecules(x) * y





def all_nucleotide(values):

    length = [2 * amplicon_length, amplicon_length, amplicon_length, primer_length, primer_length, amplicon_length + primer_length, amplicon_length + primer_length, 0, amplicon_length + primer_length, amplicon_length + primer_length, 1, extended_primer, extended_primer, amplicon_length + extended_primer,  amplicon_length + extended_primer, amplicon_length + extended_primer, amplicon_length + extended_primer]

    total_number = 0

    for i in range(len(values)):


        total_number = total_number + total_molecule_length(values[i], length[i])



    return total_number



def u_molar_concentration(values):


    length = [2 * amplicon_length, amplicon_length, amplicon_length, primer_length, primer_length, amplicon_length + primer_length, amplicon_length + primer_length, 0, amplicon_length + primer_length, amplicon_length + primer_length, 1, extended_primer, extended_primer, amplicon_length + extended_primer,  amplicon_length + extended_primer, amplicon_length + extended_primer, amplicon_length + extended_primer]

    return sum(values[i] * length[i] for i in range(len(values)))





def taq_nt_per_s(temperature):

    """This function calculates the speed of Taq DNA ploymerase, temperature is given in Kelvin"""


    taq_temperature = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(90), (90 - 0 + 1))  # x coordinates

    temperature_data = [celsius_to_Kelvin(22), celsius_to_Kelvin(37), celsius_to_Kelvin(55), celsius_to_Kelvin(70), celsius_to_Kelvin(75), celsius_to_Kelvin(80), celsius_to_Kelvin(90)]  # temperatures with known taq polymerase rate

    taq_rate_data = [0.25, 1.5, 24, 60, 150, 150, 0]           # taq polymerase rate at given temperatures

    temp_interpolar = np.interp(taq_temperature, temperature_data, taq_rate_data)        # calculating the taq polymerase rate between 0 and 90 degrees

    #plt.suptitle("Temperature-dependency of Taq polymerase", fontsize = 14)

    #plt.plot(taq_temperature, temp_interpolar)

    #plt.xlabel("Temperature (K)")

    #plt.ylabel("Incorporated nucleotide per second")

    #plt.show()


    if ((temperature >= celsius_to_Kelvin(0)) and (temperature <= celsius_to_Kelvin(90))):

        taq_index = list(taq_temperature).index(temperature)

        return temp_interpolar[taq_index]


    else:

        return 0




def denaturation(values, t, T, dGs):


    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]


    kf1 = forward_rate


    #exponent_1 = dGs[0]/(R*T)

    exponent_1 = exponent_clipping(dGs[0]/(R*T))

    kr1 = kf1 * np.exp(exponent_1)

    #print(kr1)

    kf1, kr1 = clipping(kf1, kr1)

    # if np.abs(kr1<1e-14):
    #     kr1 = 0


    #rate_den = - kr1 * S1S2 + kf1 * S1 * S2

    #rate_den = np.clip((- kr1 * S1S2 + kf1 * S1 * S2), a_min= min_clip, a_max= max_clip)

    rate_den = rate_clipping(- kr1 * S1S2 + kf1 * S1 * S2)


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


    kf2 = forward_rate



    #exponent_2a = dGs[1]/(R*T)

    #exponent_2a = np.clip((dGs[1]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_2a = exponent_clipping(dGs[1]/(R*T))

    kr2a = kf2 * np.exp(exponent_2a)

    kf2, kr2a = clipping(kf2, kr2a)

    # if np.abs(kr2a < 1e-14):
    #     kr2a = 0

    #rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2

    #rate_S1P2_bind = np.clip((kf2 * S1 * P2 - kr2a * S1P2), a_min= min_clip, a_max= max_clip)

    rate_S1P2_bind = rate_clipping(kf2 * S1 * P2 - kr2a * S1P2)

    #exponent_2b = dGs[1]/(R*T)

    #exponent_2b = np.clip((dGs[1]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_2b = exponent_clipping(dGs[1]/(R*T))

    kr2b = kf2 * np.exp(exponent_2b)

    kf2, kr2b = clipping(kf2, kr2b)

    #
    # if np.abs(kr2b < 1e-14):
    #     kr2b = 0

    #print("backwardsrate", kr2b)

    #rate_S2P1_bind = kf2 * S2 * P1 - kr2b * S2P1

    #print(rate_S2P1_bind)


    #rate_S2P1_bind = np.clip((kf2 * S2 * P1 - kr2b * S2P1), a_min= min_clip, a_max= max_clip)

    rate_S2P1_bind = rate_clipping(kf2 * S2 * P1 - kr2b * S2P1)

    #print(rate_S2P1_bind)

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


    kf5 = forward_rate

    #exponent_5a = dGs[2]/(R*T)

    #exponent_5a = np.clip((dGs[2]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_5a = exponent_clipping(dGs[2]/(R*T))

    kr5a = kf5 * np.exp(exponent_5a)

    kf5, kr5a = clipping(kf5, kr5a)

    # if np.abs(kr5a < 1e-14):
    #     kr5a = 0

    #rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    rate_S1Q2_bind = rate_clipping(kf5 * S1 * Q2 - kr5a * S1Q2)

    #exponent_5b = dGs[2]/(R*T)

    #exponent_5b = np.clip((dGs[2]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_5b = exponent_clipping(dGs[2]/(R*T))

    kr5b = kf5 * np.exp(exponent_5b)

    kf5, kr5b = clipping(kf5, kr5b)


    # if np.abs(kr5b < 1e-14):
    #     kr5b = 0

    #rate_S2Q1_bind = kf5 * S2 * Q1 - kr5b * S2Q1

    #rate_S2Q1_bind = np.clip((kf5 * S2 * Q1 - kr5b * S2Q1), a_min= min_clip, a_max= max_clip)

    rate_S2Q1_bind = rate_clipping(kf5 * S2 * Q1 - kr5b * S2Q1)

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


    kf3 = forward_rate

    #exponent_3a = dGs[3]/(R*T)

    #exponent_3a = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_3a = exponent_clipping(dGs[3]/(R*T))

    kr3a = kf3 * np.exp(exponent_3a)

    kf3, kr3a = clipping(kf3, kr3a)

    # if np.abs(kr3a<1e-14):
    #     kr3a = 0

    #rate_poly_S1P2_bind = kf3 * S1P2 * E - kr3a * S1P2E

    #rate_poly_S1P2_bind = np.clip((kf3 * S1P2 * E - kr3a * S1P2E), a_min= min_clip, a_max= max_clip)

    rate_poly_S1P2_bind = rate_clipping(kf3 * S1P2 * E - kr3a * S1P2E)

    #exponent_3b = dGs[3]/(R*T)

    #exponent_3b = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_3b = exponent_clipping(dGs[3]/(R*T))


    kr3b = kf3 * np.exp(exponent_3b)

    kf3, kr3b = clipping(kf3, kr3b)
    #
    # if np.abs(kr3b<1e-14):
    #     kr3b = 0

    #rate_poly_S2P1_bind = kf3 * S2P1 * E - kr3b * S2P1E

    #rate_poly_S2P1_bind = np.clip((kf3 * S2P1 * E - kr3b * S2P1E), a_min= min_clip, a_max= max_clip)


    rate_poly_S2P1_bind = rate_clipping(kf3 * S2P1 * E - kr3b * S2P1E)


    #enzyme_binding = - rate_poly_S1P2_bind - rate_poly_S2P1_bind

    #enzyme_binding = np.clip((- rate_poly_S1P2_bind - rate_poly_S2P1_bind), a_min= min_clip, a_max= max_clip)

    enzyme_binding = rate_clipping(- rate_poly_S1P2_bind - rate_poly_S2P1_bind)


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


    kf4 = forward_rate

    #exponent_4a = dGs[3]/(R*T)

    #exponent_4a = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_4a = exponent_clipping(dGs[3]/(R*T))


    kr4a = kf4 * np.exp(exponent_4a)

    kf4, kr4a = clipping(kf4, kr4a)

    # if np.abs(kr4a<1e-14):
    #     kr4a = 0

    #rate_poly_S1Q2_bind = kf4 * S1Q2 * E - kr4a * S1Q2E

    #rate_poly_S1Q2_bind = np.clip((kf4 * S1Q2 * E - kr4a * S1Q2E), a_min= min_clip, a_max= max_clip)

    rate_poly_S1Q2_bind = rate_clipping(kf4 * S1Q2 * E - kr4a * S1Q2E)




    #exponent_4b = dGs[3]/(R*T)

    #exponent_4b = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_4b = exponent_clipping(dGs[3]/(R*T))

    kr4b = kf4 * np.exp(exponent_4b)

    kf4, kr4b = clipping(kf4, kr4b)

    # if np.abs(kr4b<1e-14):
    #     kr4b = 0

    #rate_poly_S2Q1_bind = kf4 * S2Q1 * E - kr4b * S2Q1E

    #rate_poly_S2Q1_bind = np.clip((kf4 * S2Q1 * E - kr4b * S2Q1E), a_min= min_clip, a_max= max_clip)

    rate_poly_S2Q1_bind = rate_clipping(kf4 * S2Q1 * E - kr4b * S2Q1E)

    #enzyme = - rate_poly_S1Q2_bind - rate_poly_S2Q1_bind

    #enzyme = np.clip((- rate_poly_S1Q2_bind - rate_poly_S2Q1_bind), a_min= min_clip, a_max= max_clip)

    enzyme = rate_clipping(- rate_poly_S1Q2_bind - rate_poly_S2Q1_bind)


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

    #if (total_molecule_length(dNTP, 1)) >= (2 * n):


        #primer_ext_1.counter += 1


    #   S1P2E + n * dNTP -> S1Q2E
    #   S2P1E + n * dNTP -> S2Q1E

    #ce = 1000   # [dNTP/s] concentration of polymerase enzyme

    ce = taq_nt_per_s(T)

    #rate_ext_1 = (ce / n) * S1P2E * dNTP

    #rate_ext_1 = np.clip(((ce / n) * S1P2E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_1 = rate_clipping((ce / n) * S1P2E * dNTP)



    #rate_ext_2 = (ce / n) * S2P1E * dNTP

    #rate_ext_2 = np.clip(((ce / n) * S2P1E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_2 = rate_clipping((ce / n) * S2P1E * dNTP)

    #nucleotide = - n * rate_ext_1 - n * rate_ext_2

    #nucleotide = np.clip((- n * rate_ext_1 - n * rate_ext_2), a_min= min_clip, a_max= max_clip)

    nucleotide = rate_clipping(- n * rate_ext_1 - n * rate_ext_2)


    y = np.zeros(17)

    y[8] = - rate_ext_1  # concentration of S1P2E
    y[9] = - rate_ext_2  # concentration of S2P1E
    y[10] = nucleotide  # concentration of dNTP
    y[13] = rate_ext_1  # concentration of S1Q2E
    y[14] = rate_ext_2  # concentration of S2Q1E

    return y

    #
    # else:
    #
    #     y = np.zeros(17)
    #
    #
    #     y[8] = S1P2E    # concentration of S1P2E
    #     y[9] = S2P1E   # concentration of S2P1E
    #     y[10] = dNTP   # concentration of dNTP
    #     #y[13] = 0    # concentration of S1Q2E
    #     #y[14] = 0   # concentration of S2Q1E




def primer_ext_2(values, t, T, dGs):



    S1Q2E = values[13]
    S2Q1E = values[14]
    dNTP = values[10]

    #if total_molecule_length(dNTP, 1) >= (2 * extended_length):

       # primer_ext_2.counter += 1


        #ce_Q = 1000         # enzyme concentration

    ce_Q = taq_nt_per_s(T)

    # reaction: S1Q2E + extended_length * dNTP ---> S1S2 + E

    #rate_ext_Q1 = (ce_Q / extended_length) * S1Q2E * dNTP

    #rate_ext_Q1 = np.clip(((ce_Q / extended_length) * S1Q2E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_Q1 = rate_clipping((ce_Q / extended_length) * S1Q2E * dNTP)

    # reaction: S2Q1E + extended_length * dNTP ---> S1S2 + E

    #rate_ext_Q2 = (ce_Q / extended_length) * S2Q1E * dNTP

    #rate_ext_Q2 = np.clip(((ce_Q / extended_length) * S2Q1E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_Q2 = rate_clipping((ce_Q / extended_length) * S2Q1E * dNTP)

    #nucleotide_Q = - extended_length * rate_ext_Q1 - extended_length * rate_ext_Q2

    #nucleotide_Q = np.clip((- extended_length * rate_ext_Q1 - extended_length * rate_ext_Q2), a_min= min_clip, a_max= max_clip)

    nucleotide_Q = rate_clipping(- extended_length * rate_ext_Q1 - extended_length * rate_ext_Q2)

    # primer_ext_2.counter += 1

    y = np.zeros(17)

    y[0] = rate_ext_Q1 + rate_ext_Q2
    y[7] = rate_ext_Q1 + rate_ext_Q2
    y[10] = nucleotide_Q  # concentration of dNTP
    y[13] = - rate_ext_Q1  # concentration of S1Q2E
    y[14] = - rate_ext_Q2  # concentration of S2Q1E

    return y

    # else:
    #
    #     y = np.zeros(17)
    #
    #     y[10] = dNTP                           # concentration of dNTP
    #     y[13] = S1Q2E                         # concentration of S1Q2E
    #     y[14] = S2Q1E                          # concentration of S2Q1E
    #
    #     return y






def taq_denaturation(values, t, T, dGs):

    E = values[7]

    y = np.zeros(17)

    rate = 0.0001           #???

    #rate = 0

    taq_denaturation_T = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(110), 111)

    temperature_data = [celsius_to_Kelvin(0), celsius_to_Kelvin(90), celsius_to_Kelvin(110)]

    taq_denaturation = [0, 0, rate * E]

    temp_interpolar = np.interp(taq_denaturation_T, temperature_data, taq_denaturation)        # calculating the taq polymerase rate between 0 and 90 degrees

    # plt.title("Denaturation of Taq polymerase")
    #
    # plt.plot(taq_denaturation_T, temp_interpolar)
    #
    # plt.xlabel("Temperature (Kelvin)")
    # plt.ylabel("Rate of Taq polymerase denaturation")
    #
    # plt.show()

    if T > celsius_to_Kelvin(90):

        y[7] = - rate * E

    return y




def PCR_reaction(values, t, T, dGs):  # not using rate clipping cos the array data

    summary = denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + primer_binding_2(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + primer_ext_2(values, t, T, dGs) + taq_denaturation(values, t, T, dGs)

    summary = np.clip(summary, a_min= min_clip, a_max= max_clip)

    # summary[np.abs(summary) < 1e-14] = 0

    return summary



def S1S2_dS(values):


    functions = [denaturation, primer_binding_1, polymerase_binding_1, primer_ext_1, primer_ext_2, polymerase_binding_2, primer_binding_2]

    dGs = [0 for x in range(4)]

    #var_dS = [-52.7184e-3, -52.7184e-2, -52.7184e-1, -52.7184, -2.5 ]

    var_dS = [-5, - 3, - 2, - 1, - 0.5]

    #var_dS = [-1e-2, -1e-1, -1, -1e+1, -1.5e+1]

    time_den_S1S2 = np.linspace(0, tden, tden)

    S1S2_den_result = []


    S1S2_denature_temp = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(110), 111)     # Temperature scale between 0 and 110 degree Celsius

    var_dS_result = np.zeros((len(var_dS), len(S1S2_denature_temp)))


    for x in range(len(var_dS)):


        dS = var_dS[x]



        for i in range(len(S1S2_denature_temp)):


            dGs[0] = (Tm_S1S2 - S1S2_denature_temp[i]) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+14)

            integration_den = odeint(functions[0], values, time_den_S1S2, args=(S1S2_denature_temp[i], dGs))

            #S1S2_den_result.append(integration_den[-1, 0])       # storing the concentration of S1S2 at the end of denaturation in this variable

            #print(len(S1S2_den_result))

            var_dS_result[x, i] = integration_den[-1, 0]

            values = initial_fixed                               # the initial concentrations will be used in the next cycle at a different temperature






    # half_value = initial_fixed[0]/2
    #
    # closest_value = 100             # initialized the smallest difference in a high number
    #
    # for i in range(len(S1S2_den_result)):
    #
    #     if closest_value > np.abs(S1S2_den_result[i] - half_value):
    #
    #         closest_value = S1S2_den_result[i]
    #
    #
    # index_temp = S1S2_den_result.index(closest_value)
    #
    # print(" Actual melting temperature:", S1S2_denature_temp[index_temp],"\n","Half-concentration of the initial S1S2:", initial_fixed[0]/2,"\n","The closest S1S2 concentration based on integer temperature (Celsius):", closest_value)

    plt.title("S1S2 concentration with different dS values\nOn the temperature scale from 273.15 to 383.15 K (0 to 110 °C)", FontSize= 16, FontWeight = "bold", position=(0.5, 1.05))


    for x in range(len(var_dS)):


        plt.plot(S1S2_denature_temp, var_dS_result[x, ])



    plt.axvline(x=Tm_S1S2, color = "b", linestyle= ":")

    #plt.legend(["Tm"])

    plt.axhline(y=(initial_fixed[0]/2), color = "k", linestyle= ":")

    plt.xlabel("Temperature (K)", FontSize= 13, FontWeight = "bold", position=(0.9,-1))


    plt.ylabel("Concentration of S1S2 (umol)", FontSize= 13, FontWeight = "bold", position=(0,0.6))

    #curve_legend = ["dS =", var_dS[x]]

    #plt.legend(["S1S2", "Tm_S1S2", "Half-concentration"])

    plt.legend(["dS = " + str(var_dS[0]), "dS = " + str(var_dS[1]), "dS = " + str(var_dS[2]), "dS = " + str(var_dS[3]), "dS = " + str(var_dS[4]), "Tm_S1S2", "Half-concentration of S1S2"], loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    plt.tight_layout()

    plt.show()


    #return S1S2_denature_temp[index_temp]


def dS1S2_dS_Tden(values):


    functions = [denaturation, primer_binding_1, polymerase_binding_1, primer_ext_1, primer_ext_2, polymerase_binding_2, primer_binding_2]

    dGs = [0 for x in range(4)]

    #var_dS = [-1, -2, -3, -4, -5]

    var_dS = [-5, - 3, - 2, - 1, - 0.5]

    #var_dS = [- 2, - 0.99, - 0.5, - 0.3, - 0.1]

    var_S1S2 = np.linspace(1e-8, 1e+3, 111)

    time_den_S1S2 = np.linspace(0, tden, tden)
    var_dS_result = np.zeros((len(var_dS), len(var_S1S2)))


    for x in range(len(var_dS)):


        dS = var_dS[x]



        for i in range(len(var_S1S2)):

            values[0] = var_S1S2[i]


            dGs[0] = (Tm_S1S2 - Tden) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+14)

            integration_den = odeint(functions[0], values, time_den_S1S2, args=(Tden, dGs))

            var_dS_result[x, i] = integration_den[-1, 0]

            values = initial_fixed                               # the initial concentrations will be used in the next cycle at a different temperature




    #print(var_dS_result)


    plt.title("S1S2 concentration with different dS values\nThe initial S1S2 concentration varies between 1e-8 and 1e+3 Micromolar", FontSize= 14, FontWeight = "bold")


    for x in range(len(var_dS)):


        plt.plot(var_S1S2, (var_dS_result[x, ]/var_S1S2)*100)


    plt.xlabel("Initial concentration of S1S2", FontSize= 11, FontWeight = "bold")

    plt.ylabel("Percentage of S1S2 concentration (umol) \n at the end of denaturation \n compared to the initial S1S2 concentration", FontSize= 11, FontWeight = "bold")


    plt.legend(["dS = " + str(var_dS[0]), "dS = " + str(var_dS[1]), "dS = " + str(var_dS[2]), "dS = " + str(var_dS[3]), "dS = " + str(var_dS[4])], loc='upper left', prop={'size':10}, bbox_to_anchor=(1,1))

    plt.tight_layout()

    plt.show()


def individual_integration(values):


    functions = [denaturation, primer_binding_1, polymerase_binding_1, primer_ext_1, polymerase_binding_2, primer_binding_2, primer_ext_2]


    functions_name = ["Denaturation", "Primer binding 1", "Polymerase binding 1", "Primer extension 1", "Polymerase binding 2", "Primer binding 2", "Primer extension 2"]

    concentration = np.empty((number_time_points, 17))

    dGs = [0 for x in range(4)]

    difference_umol = []


    for n in range(len(functions)):



        before_nt_number = all_nucleotide(values)

        all_umol_1 = u_molar_concentration(values)


        for i in range(number_cycles):


            dGs[0] = (Tm_S1S2 - Tden) * dS

            dGs[1] = (Tm_primer - Tden) * dS

            dGs[2] = (Tm_extended_primer - Tden) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tden * dS

            dGs[3] = (Tm_enzyme - Tden) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_den", dGs)


            integration_den = odeint(functions[n], values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs))

            concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den


            dGs[0] = (Tm_S1S2 - Tanneal) * dS

            dGs[1] = (Tm_primer - Tanneal) * dS

            dGs[2] = (Tm_extended_primer - Tanneal) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tanneal * dS

            dGs[3] = (Tm_enzyme - Tanneal) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_anneals", dGs)



            integration_anneal = odeint(functions[n], integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs))

            concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal



            dGs[0] = (Tm_S1S2 - Text) * dS

            dGs[1] = (Tm_primer - Text) * dS

            dGs[2] = (Tm_extended_primer - Text) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Text * dS

            dGs[3] = (Tm_enzyme - Text) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_text", dGs)

            integration_ext = odeint(functions[n], integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs))

            concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


            values = integration_ext[-1]

        print("The concentration of the 17 species after", functions_name[n], "is added to the process individually:", values)


        after_nt_number = all_nucleotide(values)

        all_umol_2 = u_molar_concentration(values)

        difference_umol.append(all_umol_1 - all_umol_2)

        print("The difference in nt number after", functions_name[n], ":", before_nt_number - after_nt_number, "\n")
        print("The difference in micromol after", functions_name[n], ":", all_umol_1 - all_umol_2, "\n")


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

    #print(difference_umol)

    plt.figure(2)

    plt.title("Difference of micromolar concentration after each function is added individually", FontSize= 18, FontWeight = "bold", position=(0.5, 1.05))

    plt.plot(range(len(functions)), difference_umol)

    plt.xticks(range(len(functions)), functions_name)

    plt.xlabel("Name of the individually added function", FontSize= 13, FontWeight = "bold", position=(0.9, 0))
    plt.ylabel("Difference in micromolar concentration", FontSize= 13, FontWeight = "bold", position=(0,0.8))

    plt.show()


    #the created values will be used in the next function


    #print("The concentration of the 17 species at the end of the process:", values)

    return values





def a(values, t, T, dGs):

    return denaturation(values, t, T, dGs)


def b(values, t, T, dGs):

    return denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs)


def c(values, t, T, dGs):

    return denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs)


def d(values, t, T, dGs):

    return denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + primer_ext_1(values, t, T, dGs)



def e(values, t, T, dGs):

    return denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs)




def f(values, t, T, dGs):

    return denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) + primer_binding_2(values, t, T, dGs)




def g(values, t, T, dGs):

    return denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) +  primer_binding_2(values, t, T, dGs) + primer_ext_2(values, t, T, dGs)





def iterative_integration(values, number):



    all_summaries = [a,b,c, d,e,f,g]

    difference_umol = []

    for n in range(number):


        concentration = np.empty((number_time_points, 17))

        dGs = [0 for x in range(4)]

        primer_ext_1.counter = 0

        primer_ext_2.counter = 0

        before_nt_number = all_nucleotide(values)

        all_umol_1 = u_molar_concentration(values)


        for i in range(number_cycles):


            dGs[0] = (Tm_S1S2 - Tden) * dS

            dGs[1] = (Tm_primer - Tden) * dS

            dGs[2] = (Tm_extended_primer - Tden) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tden * dS

            dGs[3] = (Tm_enzyme - Tden) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_den", dGs)


            integration_den = odeint(all_summaries[n], values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs))

            concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den

            dGs[0] = (Tm_S1S2 - Tanneal) * dS

            dGs[1] = (Tm_primer - Tanneal) * dS

            dGs[2] = (Tm_extended_primer - Tanneal) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tanneal * dS

            dGs[3] = (Tm_enzyme - Tanneal) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_anneals", dGs)

            integration_anneal = odeint(all_summaries[n], integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs))

            concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal

            dGs[0] = (Tm_S1S2 - Text) * dS

            dGs[1] = (Tm_primer - Text) * dS

            dGs[2] = (Tm_extended_primer - Text) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Text * dS

            dGs[3] = (Tm_enzyme - Text) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_text", dGs)

            integration_ext = odeint(all_summaries[n], integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs))

            concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


            values = integration_ext[-1]

        print("The concentration of the 17 species after", functions_name[n],  "is added to the process iteratively:", values)


        after_nt_number = all_nucleotide(values)

        all_umol_2 = u_molar_concentration(values)

        difference_umol.append(all_umol_1 - all_umol_2)

        print("The difference in nt number after", functions_name[n], ":", before_nt_number - after_nt_number, "\n")

        print("The difference in micromol after", functions_name[n], ":", all_umol_1 - all_umol_2, "\n")


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

        values = initial_fixed


    plt.figure(2)

    plt.title("Difference of micromolar concentration after each function is added iteratively",  FontSize= 18, FontWeight = "bold", position=(0.5, 1.05))

    plt.plot(range(number), difference_umol)

    plt.xticks(range(number), functions_name)

    plt.xlabel("Name of the iteratively added function", FontSize= 13, FontWeight = "bold", position=(0.9, 0))
    plt.ylabel("Difference in micromolar concentration", FontSize= 13, FontWeight = "bold", position=(0,0.8))

    plt.show()


    return integration_ext[-1]



def only_one_integration(values, number):

    number = number - 1


    functions_plus = [denaturation, primer_binding_1, polymerase_binding_1, primer_ext_1, polymerase_binding_2, primer_binding_2, primer_ext_2, PCR_reaction]

    functions_plus_name = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer_ext_2", "PCR_reaction"]


    concentration = np.empty((number_time_points, 17))

    dGs = [0 for x in range(4)]



    before_nt_number = all_nucleotide(values)

    all_umol_1 = u_molar_concentration(values)


    for i in range(number_cycles):


        dGs[0] = (Tm_S1S2 - Tden) * dS

        dGs[1] = (Tm_primer - Tden) * dS

        dGs[2] = (Tm_extended_primer - Tden) * dS

        #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tden * dS

        dGs[3] = (Tm_enzyme - Tden) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        #print("dGs_den", dGs)




        integration_den = odeint(functions_plus[number], values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs))

        concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den

        dGs[0] = (Tm_S1S2 - Tanneal) * dS

        dGs[1] = (Tm_primer - Tanneal) * dS

        dGs[2] = (Tm_extended_primer - Tanneal) * dS

        #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tanneal * dS

        dGs[3] = (Tm_enzyme - Tanneal) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        #print("dGs_anneals", dGs)

        integration_anneal = odeint(functions_plus[number], integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs))

        concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal

        dGs[0] = (Tm_S1S2 - Text) * dS

        dGs[1] = (Tm_primer - Text) * dS

        dGs[2] = (Tm_extended_primer - Text) * dS

        #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Text * dS

        dGs[3] = (Tm_enzyme - Text) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        #print("dGs_text", dGs)

        integration_ext = odeint(functions_plus[number], integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs))

        concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


        values = integration_ext[-1]




    print("The concentration of the 17 species at the end of", functions_plus_name[number], ":", values)

    after_nt_number = all_nucleotide(values)

    all_umol_2 = u_molar_concentration(values)

    print("The difference in nt number after", functions_plus_name[number], ":", before_nt_number - after_nt_number, "\n")

    print("The difference in micromolar concentration after", functions_plus_name[number], ":", all_umol_1 - all_umol_2, "\n")


    #Plotting the concentrations over time

    plt.figure(1)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)

    plots = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(10):

        plt.subplot(2, 5, i+1)


        plt.plot(time, concentration[:, plots[i]])

        #plt.ylim([0, y_top_limit[i]])


        plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time")
        plt.ylabel("Concentration")




    plt.show()




    return values




if __name__ == '__main__':


    #print(u_molar_concentration(values))



    #S1S2_dS(values)

    #dS1S2_dS_Tden(values)

    #print(taq_denaturation(values, Tanneal))




    #individual_integration(values)

    #iterative_integration(values, 7)






    only_one_integration(values, 8)




    #print(values)

    #print(u_molar_concentration(values))









