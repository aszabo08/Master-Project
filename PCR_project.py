
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


new_species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1", "S1M2", "S1M2E", "S1N2E", "S1L2", "L2", "S1N2", "L2P1", "L2P1E", "L2Q1E", "L1L2", "L2Q1", "L1", "L1P2", "L1P2E", "L1Q2E", "L1Q2"]


values = [0 for i in range(33)]


#values[0] = 3.0769e-2       # concentration of plasmid (S1S2) in uM

values[0] = 0.000445632798573975

values[3] = 0.5       # concentration of P1 in uM
values[4] = 0.5

values[7] = 0.02

values[10] = 200


initial_dNTP = values[10]

initial_fixed = values

functions_name = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer_ext_2"]


T_initial_den = 371.15

Tden = 371.15
Tanneal = 330.15
Text = 340.15
T_cooling_down = 298.15             # 25 degree celsius


t_initial_den = 30

tden = 10                   # seconds
tanneal = 0                # seconds
text = 20                   # seconds
t_cooling_down = 300


total = tden + tanneal + text
number_cycles = 32
steps = 1


# amplicon_length = 1000
# primer_length = 15
# n = 10

amplicon_length = 340
primer_length = 28
n = 10


extended_primer = primer_length + n
extended_length = amplicon_length - extended_primer

time = np.linspace(0, t_initial_den + number_cycles * (tden + tanneal + text) + t_cooling_down, number_cycles * (tden + tanneal + text) * steps + (t_cooling_down + t_initial_den) * steps)  # for every second "steps" points are distinguished

number_time_points = total * steps * number_cycles + (t_initial_den + t_cooling_down) * steps

enzyme_type = 'q5'



# original Tm -es

#
# Tm_S1S2 = 363.15                # 90 degrees
#
# Tm_primer = 308.15              # 35 degrees        # from neb tm calculator: 15nt, 20 % GC content, 35 celius TM: AATTTAACGAATTCA
#
# Tm_extended_primer = 313.15     # 40 degree
#
# Tm_enzyme = 353.15              # 80 degree

dS = - 2


max_exponent = 15    # 13
min_clip = -1e+15       #18
max_clip = 1e+15

forward_rate = 1


R = 8.314e-3        # Gas contant in  J / K mol

#Tm_primer = 320.15      # 47 deggree

#Tm_extended_primer = 337.15     #64 degree


# Tm_primer = 343.15
#
# Tm_extended_primer = 349.15


Tm_primer = 345.15
#
Tm_extended_primer = 348.15



Tmax = 373.15           # 100 degree


K = (primer_length * extended_primer * (Tm_extended_primer - Tm_primer)) / (extended_primer * Tm_primer - primer_length * Tm_extended_primer)


dH = (Tm_primer * dS * (primer_length + K)) / (Tmax * primer_length)


dH_check = (Tm_extended_primer * dS * ( extended_primer + K)) / ( Tmax * extended_primer)

Tm_S1S2 = (Tmax * amplicon_length * dH) / ((amplicon_length + K) * dS)
#
Tm_enzyme1 = (Tmax * 1 * dH) / ((1 + K) * dS)

Tm_enzyme2 = dH/dS

Tm_enzyme = 356.15             # 85 degree  not calculated!










def plot_attributes():


    plt.xlabel("Time (s)")
    plt.ylabel("Concentration (uM)")
    matplotlib.pyplot.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace=0.2)




def celsius_to_Kelvin(x):

    return x + 273.15


def clipping(kf, kr):

    return (kf / kr) * np.clip(kr, a_min= min_clip, a_max= max_clip), np.clip(kr, a_min= min_clip, a_max= max_clip)


def rate_clipping(x):

    #if x != np.clip(x, a_min= min_clip, a_max= max_clip):

        #print("Clipped rate value:", x)


    return np.clip(x, a_min= min_clip, a_max= max_clip)


def exponent_clipping(x):

    #if x != np.clip(x, a_min= None, a_max= max_exponent):

        #print("Clipped exponent value:", x)


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


    length = [2 * amplicon_length, amplicon_length, amplicon_length, primer_length, primer_length, amplicon_length + primer_length, amplicon_length + primer_length, 0, amplicon_length + primer_length, amplicon_length + primer_length, 1, extended_primer, extended_primer, amplicon_length + extended_primer,  amplicon_length + extended_primer, amplicon_length + extended_primer, amplicon_length + extended_primer, ]

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


    if ((temperature >= celsius_to_Kelvin(22)) and (temperature <= celsius_to_Kelvin(90))):

        taq_index = list(taq_temperature).index(temperature)

        return temp_interpolar[taq_index]


    else:

        return 0





def polymerase_nt_per_s(temperature):

    """This function calculates the speed of DNA ploymerase, temperature is given in Kelvin"""


    #print(enzyme_type)

    temperature_scale = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(95), (95 - 0 + 1))  # x coordinates


    taq_temperature_data = [celsius_to_Kelvin(20), celsius_to_Kelvin(22), celsius_to_Kelvin(37), celsius_to_Kelvin(55), celsius_to_Kelvin(70), celsius_to_Kelvin(75), celsius_to_Kelvin(80), celsius_to_Kelvin(90)]  # temperatures with known taq polymerase rate

    taq_rate_data = [0, 0.25, 1.5, 24, 60, 150, 150, 0]           # taq polymerase rate at given temperatures

    taq_temp_interpolar = np.interp(temperature_scale, taq_temperature_data, taq_rate_data)        # calculating the taq polymerase rate between 0 and 90 degrees


    q5_temperature_data = [celsius_to_Kelvin(25), celsius_to_Kelvin(30), celsius_to_Kelvin(45), celsius_to_Kelvin(60), celsius_to_Kelvin(70), celsius_to_Kelvin(75), celsius_to_Kelvin(80), celsius_to_Kelvin(85), celsius_to_Kelvin(90), celsius_to_Kelvin(93) ]  # temperatures with known taq polymerase rate

    q5_rate_data = [0, 4.63, 9.26, 46.30, 83.33, 106.48, 85.19, 78.70, 13.89, 0]           # taq polymerase rate at given temperatures

    q5_temp_interpolar = np.interp(temperature_scale, q5_temperature_data, q5_rate_data)        # calculating the taq polymerase rate between 0 and 90 degrees





    #
    # plt.suptitle("Temperature-dependency of Taq polymerase", fontsize = 14)
    #
    # plt.plot(temperature_scale, taq_temp_interpolar)
    #
    # plt.plot(temperature_scale, q5_temp_interpolar)
    #
    # plt.xlabel("Temperature (K)")
    #
    # plt.ylabel("Incorporated nucleotide per second")
    #
    # plt.show()



    if enzyme_type == 'taq':


        if ((temperature >= celsius_to_Kelvin(20)) and (temperature <= celsius_to_Kelvin(90))):

            taq_index = list(temperature_scale).index(temperature)

            return taq_temp_interpolar[taq_index]


        else:

            return 0


    elif enzyme_type == 'q5':


        if ((temperature >= celsius_to_Kelvin(25)) and (temperature <= celsius_to_Kelvin(93))):

            q5_index = list(temperature_scale).index(temperature)

            return q5_temp_interpolar[q5_index]


        else:

            return 0










def denaturation(values, t, T, dGs):


    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]


    kf1 = forward_rate


    exponent_1 = exponent_clipping(dGs[0]/(R*T))

    kr1 = kf1 * np.exp(exponent_1)

    kf1, kr1 = clipping(kf1, kr1)

    rate_den = rate_clipping(- kr1 * S1S2 + kf1 * S1 * S2)


    y = np.zeros(33)

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

    y = np.zeros(33)


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

    y = np.zeros(33)


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


    y = np.zeros(33)

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


    y = np.zeros(33)

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

    ce = polymerase_nt_per_s(T)

    #rate_ext_1 = (ce / n) * S1P2E * dNTP

    #rate_ext_1 = np.clip(((ce / n) * S1P2E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_1 = rate_clipping((ce / n) * S1P2E * dNTP)



    #rate_ext_2 = (ce / n) * S2P1E * dNTP

    #rate_ext_2 = np.clip(((ce / n) * S2P1E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_2 = rate_clipping((ce / n) * S2P1E * dNTP)

    #nucleotide = - n * rate_ext_1 - n * rate_ext_2

    #nucleotide = np.clip((- n * rate_ext_1 - n * rate_ext_2), a_min= min_clip, a_max= max_clip)

    nucleotide = rate_clipping(- n * rate_ext_1 - n * rate_ext_2)


    y = np.zeros(33)

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

    ce_Q = polymerase_nt_per_s(T)

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


    product = rate_clipping(rate_ext_Q1 + rate_ext_Q2)




    # primer_ext_2.counter += 1

    y = np.zeros(33)


    #clipping for sums? --- product


    y[0] = product
    y[7] = product

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

    y = np.zeros(33)

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

    #summary = rate_clipping(summary)

    # summary[np.abs(summary) < 1e-14] = 0

    return summary




def enzyme_dissociation(values, t, T, dGs):

    return primer_binding_1(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs)





low_concentration = [0 for i in range(33)]

low_concentration[0] = 45

low_concentration[5], low_concentration[6] = 8, 8

low_concentration[8], low_concentration[9] = 0.0025, 0.0025

low_concentration[15], low_concentration[16] = 1.5e-8, 1.5e-8


high_concentration = [0 for i in range(33)]

high_concentration[0] = 4000

high_concentration[5], high_concentration[6] = 60, 60

high_concentration[8], high_concentration[9] = 0.0005, 0.0005

high_concentration[15], high_concentration[16] = 2.5e-8, 2.5e-8

overall_concentration = [low_concentration, high_concentration]

initial_overall = overall_concentration






def dS_change_4_species2(overall_concentration):





    functions = [denaturation, primer_binding_1, primer_binding_2, enzyme_dissociation]

    dGs = [0 for x in range(4)]

    #var_dS = [-52.7184e-3, -52.7184e-2, -52.7184e-1, -52.7184, -2.5 ]

    var_dS = [-10, -5, -1, -0.5]

    #var_dS = [-5, - 3, - 2, - 1, - 0.5]

    #var_dS = [-1e-2, -1e-1, -1, -1e+1, -1.5e+1]



    time_all = np.linspace(0, tden + tanneal + text, tden + tanneal + text)


    temperature_scale = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(110), 111)     # Temperature scale between 0 and 110 degree Celsius

    var_dS_S1S2_result = np.zeros((len(overall_concentration), len(var_dS), len(temperature_scale)))

    var_dS_primer_result = np.zeros((len(overall_concentration), len(var_dS), len(temperature_scale)))

    var_dS_ext_primer_result = np.zeros((len(overall_concentration), len(var_dS), len(temperature_scale)))

    var_dS_enzyme_result = np.zeros((len(overall_concentration), len(var_dS), len(temperature_scale)))


    for e in range(len(overall_concentration)):



        #values = overall_concentration[e]



        for x in range(len(var_dS)):


            dS = var_dS[x]



            for i in range(len(temperature_scale)):



                dGs[0] = (Tm_S1S2 - temperature_scale[i]) * dS

                dGs[1] = (Tm_primer - temperature_scale[i]) * dS

                dGs[2] = (Tm_extended_primer - temperature_scale[i]) * dS

                dGs[3] = (Tm_enzyme - temperature_scale[i]) * dS


                dGs = np.clip(dGs, a_min=None, a_max=1e+14)


                integration_den = odeint(functions[0], overall_concentration[e], time_all, args=(temperature_scale[i], dGs),  mxstep=5000000)

                var_dS_S1S2_result[e, x, i] = integration_den[-1, 0]

                overall_concentration = initial_overall



                integration_primer = odeint(functions[1], overall_concentration[e], time_all, args=(temperature_scale[i], dGs),  mxstep=5000000)

                var_dS_primer_result[e, x, i] = integration_primer[-1, 5]

                overall_concentration = initial_overall



                integration_ext_primer = odeint(functions[2], overall_concentration[e], time_all, args=(temperature_scale[i], dGs),  mxstep=5000000)

                var_dS_ext_primer_result[e, x, i] = integration_ext_primer[-1, 15]

                overall_concentration = initial_overall                              # the initial concentrations will be used in the next cycle at a different temperature




                integration_enzyme = odeint(functions[3], overall_concentration[e], time_all, args=(temperature_scale[i], dGs),  mxstep=5000000)

                var_dS_enzyme_result[e, x, i] = integration_enzyme[-1, 8]

                overall_concentration = initial_overall







    plt.figure(1)

    plt.title("Low (-) and high (:) S1S2 initial concentration with different dS values", FontSize= 16, FontWeight = "bold", position=(0.5, 1.05))

    style_curve = ['-',':']

    colour_curve = ['C1', 'C2', 'C3','C4']


    for i in range(len(overall_concentration)):


        for x in range(len(var_dS)):


            plt.plot(temperature_scale, var_dS_S1S2_result[i, x, ]/initial_overall[i][0]*100, style_curve[i], color = colour_curve[x], label = "dS = " + str(var_dS[x]))     #/initial_overall[e][0]*100



    plt.axvline(x=Tm_S1S2, color = 'black', linestyle= "--", label = "Tm_S1S2")

    #plt.legend(["Tm"])

    #plt.axhline(y=(initial_fixed[0]/2), color = "k", linestyle= ":")

    plt.xlim(320, None)

    plt.xlabel("Temperature (K)", FontSize= 13, FontWeight = "bold", position=(0.9,-1))


    #plt.ylabel("Concentration of S1S2 (uM)", FontSize= 13, FontWeight = "bold", position=(0,0.6))

    plt.ylabel("Percentage of S1S2 concentration (umol) \n at the end of denaturation \n compared to the initial S1S2 concentration", FontSize= 11, FontWeight = "bold")

    #curve_legend = ["dS =", var_dS[x]]

    #plt.legend(["S1S2", "Tm_S1S2", "Half-concentration"])

    #plt.legend(["dS = " + str(var_dS[0]), "dS = " + str(var_dS[1]), "dS = " + str(var_dS[2]), "dS = " + str(var_dS[3]), "dS = " + str(var_dS[4]), "Tm_S1S2", "Half-concentration of S1S2"], loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    plt.legend(loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    #plt.tight_layout()



    plt.figure(2)

    plt.title("Low (-) and high (:) initial primer concentration with different dS values", FontSize= 16, FontWeight = "bold", position=(0.5, 1.05))


    for i in range(len(overall_concentration)):


        for x in range(len(var_dS)):


            plt.plot(temperature_scale, var_dS_primer_result[i, x, ]/initial_overall[i][5]*100, style_curve[i], color = colour_curve[x], label = "dS = " + str(var_dS[x]))          #/initial_overall[e][5]*100



    plt.axvline(x=Tm_primer, color = 'black', linestyle= "--", label = "Tm_primer")

    #plt.legend(["Tm"])

    #plt.axhline(y=(initial_fixed[3]/2), color = "k", linestyle= ":")

    plt.xlim(None ,360)

    plt.xlabel("Temperature (K)", FontSize= 13, FontWeight = "bold", position=(0.9,-1))


    #plt.ylabel("Concentration of primer (uM)", FontSize= 13, FontWeight = "bold", position=(0,0.6))

    plt.ylabel("Percentage of primer concentration (umol) \n at the end of primer binding 1 reaction\n compared to the initial primer concentration", FontSize= 11, FontWeight = "bold")

    #curve_legend = ["dS =", var_dS[x]]

    #plt.legend(["S1S2", "Tm_S1S2", "Half-concentration"])

    #plt.legend(["dS = " + str(var_dS[0]), "dS = " + str(var_dS[1]), "dS = " + str(var_dS[2]), "dS = " + str(var_dS[3]), "dS = " + str(var_dS[4]), "Tm_primer", "Half-concentration of primer"], loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    plt.legend(loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    #plt.tight_layout()




    plt.figure(3)

    plt.title("Low (-) and high (:) initial extended primer concentration with different dS values", FontSize= 16, FontWeight = "bold", position=(0.5, 1.05))


    for i in range(len(overall_concentration)):


        for x in range(len(var_dS)):


            #rint("just var",var_dS_ext_primer_result[i, x, ])

            #print("initial", initial_overall[i][15]*100)


            plt.plot(temperature_scale, var_dS_ext_primer_result[i, x, ]/initial_overall[i][15]*100, style_curve[i], color = colour_curve[x], label = "dS = " + str(var_dS[x]))      # /initial_overall[e][15]*100



    plt.axvline(x=Tm_extended_primer, color = 'black', linestyle= "--", label = "Tm_extended_primer")

    #plt.legend(["Tm"])

    #plt.axhline(y=(initial_fixed[0]/2), color = "k", linestyle= ":")

    plt.xlabel("Temperature (K)", FontSize= 13, FontWeight = "bold", position=(0.9,-1))


    #plt.ylabel("Concentration of extended primer (uM)", FontSize= 13, FontWeight = "bold", position=(0,0.6))

    plt.ylabel("Percentage of extended primer concentration (umol) \n at the end of primer binding 2 reaction\n compared to the initial extended primer concentration", FontSize= 11, FontWeight = "bold")

    #curve_legend = ["dS =", var_dS[x]]

    #plt.legend(["S1S2", "Tm_S1S2", "Half-concentration"])

    #plt.legend(["dS = " + str(var_dS[0]), "dS = " + str(var_dS[1]), "dS = " + str(var_dS[2]), "dS = " + str(var_dS[3]), "dS = " + str(var_dS[4]), "Tm_extended_primer"], loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    plt.legend(loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    #plt.tight_layout()



    plt.figure(4)

    plt.title("Low (-) and high (:) initial enzyme concentration with different dS values", FontSize= 16, FontWeight = "bold", position=(0.5, 1.05))


    for i in range(len(overall_concentration)):


        for x in range(len(var_dS)):


            plt.plot(temperature_scale, var_dS_enzyme_result[i, x, ]/initial_overall[i][8]*100, style_curve[i], color = colour_curve[x], label = "dS = " + str(var_dS[x]))           # /initial_overall[e][8]*100



    plt.axvline(x=Tm_enzyme, color = 'black', linestyle= "--", label = "Tm_enzyme")

    #plt.legend(["Tm"])

    #plt.axhline(y=(initial_fixed[7]/2), color = "k", linestyle= ":")

    plt.xlim(320, 380)

    plt.xlabel("Temperature (K)", FontSize= 13, FontWeight = "bold", position=(0.9,-1))


    #plt.ylabel("Concentration of enzyme (uM)", FontSize= 13, FontWeight = "bold", position=(0,0.6))

    plt.ylabel("Percentage of enzyme concentration (umol) \n at the end of primer binding 1 and polymerase binding 1 reactions\n compared to the initial enzyme concentration", FontSize= 11, FontWeight = "bold")

    #curve_legend = ["dS =", var_dS[x]]

    #plt.legend(["S1S2", "Tm_S1S2", "Half-concentration"])

    #plt.legend(["dS = " + str(var_dS[0]), "dS = " + str(var_dS[1]), "dS = " + str(var_dS[2]), "dS = " + str(var_dS[3]), "dS = " + str(var_dS[4]), "Tm_enzyme", "Half-concentration of enzyme"], loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    plt.legend(loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    #plt.tight_layout()




    plt.show()


    #return S1S2_denature_temp[index_temp]






def individual_integration(values):


    functions = [denaturation, primer_binding_1, polymerase_binding_1, primer_ext_1, polymerase_binding_2, primer_binding_2, primer_ext_2]


    functions_name = ["Denaturation", "Primer binding 1", "Polymerase binding 1", "Primer extension 1", "Polymerase binding 2", "Primer binding 2", "Primer extension 2"]

    concentration = np.empty((number_time_points, 33))

    dGs = [0 for x in range(4)]

    difference_umol = []


    for n in range(len(functions)):



        #before_nt_number = all_nucleotide(values)

        #all_umol_1 = u_molar_concentration(values)


        for i in range(number_cycles):


            dGs[0] = (Tm_S1S2 - Tden) * dS

            dGs[1] = (Tm_primer - Tden) * dS

            dGs[2] = (Tm_extended_primer - Tden) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tden * dS

            dGs[3] = (Tm_enzyme - Tden) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_den", dGs)


            integration_den = odeint(functions[n], values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

            concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den


            dGs[0] = (Tm_S1S2 - Tanneal) * dS

            dGs[1] = (Tm_primer - Tanneal) * dS

            dGs[2] = (Tm_extended_primer - Tanneal) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tanneal * dS

            dGs[3] = (Tm_enzyme - Tanneal) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_anneals", dGs)



            integration_anneal = odeint(functions[n], integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)

            concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal



            dGs[0] = (Tm_S1S2 - Text) * dS

            dGs[1] = (Tm_primer - Text) * dS

            dGs[2] = (Tm_extended_primer - Text) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Text * dS

            dGs[3] = (Tm_enzyme - Text) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_text", dGs)

            integration_ext = odeint(functions[n], integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)

            concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


            values = integration_ext[-1]

        print("The concentration of the 17 species after", functions_name[n], "is added to the process individually:", values)


        #after_nt_number = all_nucleotide(values)

        #all_umol_2 = u_molar_concentration(values)

        #difference_umol.append(all_umol_1 - all_umol_2)

        #print("The difference in nt number after", functions_name[n], ":", before_nt_number - after_nt_number, "\n")
        #print("The difference in micromol after", functions_name[n], ":", all_umol_1 - all_umol_2, "\n")


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

    #plt.figure(2)

    #plt.title("Difference of micromolar concentration after each function is added individually", FontSize= 18, FontWeight = "bold", position=(0.5, 1.05))

    #plt.plot(range(len(functions)), difference_umol)

    #plt.xticks(range(len(functions)), functions_name)

    #plt.xlabel("Name of the individually added function", FontSize= 13, FontWeight = "bold", position=(0.9, 0))
    #plt.ylabel("Difference in micromolar concentration", FontSize= 13, FontWeight = "bold", position=(0,0.8))

    #plt.show()


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

    #difference_umol = []

    for n in range(number):


        concentration = np.empty((number_time_points, 33))

        dGs = [0 for x in range(4)]



        #before_nt_number = all_nucleotide(values)

        #all_umol_1 = u_molar_concentration(values)


        for i in range(number_cycles):


            dGs[0] = (Tm_S1S2 - Tden) * dS

            dGs[1] = (Tm_primer - Tden) * dS

            dGs[2] = (Tm_extended_primer - Tden) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tden * dS

            dGs[3] = (Tm_enzyme - Tden) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_den", dGs)


            integration_den = odeint(all_summaries[n], values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

            concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den

            dGs[0] = (Tm_S1S2 - Tanneal) * dS

            dGs[1] = (Tm_primer - Tanneal) * dS

            dGs[2] = (Tm_extended_primer - Tanneal) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Tanneal * dS

            dGs[3] = (Tm_enzyme - Tanneal) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_anneals", dGs)

            integration_anneal = odeint(all_summaries[n], integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)

            concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal

            dGs[0] = (Tm_S1S2 - Text) * dS

            dGs[1] = (Tm_primer - Text) * dS

            dGs[2] = (Tm_extended_primer - Text) * dS

            #dGs[3] = ((20 * Tm_primer * dS) / primer_length) - Text * dS

            dGs[3] = (Tm_enzyme - Text) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            #print("dGs_text", dGs)

            integration_ext = odeint(all_summaries[n], integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)

            concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


            values = integration_ext[-1]

        print("The concentration of the 17 species after", functions_name[n],  "is added to the process iteratively:", values)


        #after_nt_number = all_nucleotide(values)

        #all_umol_2 = u_molar_concentration(values)

        #difference_umol.append(all_umol_1 - all_umol_2)

        #print("The difference in nt number after", functions_name[n], ":", before_nt_number - after_nt_number, "\n")

        #print("The difference in micromol after", functions_name[n], ":", all_umol_1 - all_umol_2, "\n")


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


    #plt.figure(2)

    #plt.title("Difference of micromolar concentration after each function is added iteratively",  FontSize= 18, FontWeight = "bold", position=(0.5, 1.05))

    #plt.plot(range(number), difference_umol)

    #plt.xticks(range(number), functions_name)

    #plt.xlabel("Name of the iteratively added function", FontSize= 13, FontWeight = "bold", position=(0.9, 0))
    #plt.ylabel("Difference in micromolar concentration", FontSize= 13, FontWeight = "bold", position=(0,0.8))

    #plt.show()


    return integration_ext[-1]



def PCR_integration(values):

    concentration = np.empty((number_time_points, 33))          # len new species!

    dGs = [0 for x in range(4)]

    #before_nt_number = all_nucleotide(values)

    #all_umol_1 = u_molar_concentration(values)







    for i in range(number_cycles):

        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tden * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tden * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tden * dS)

        dGs[3] = (Tm_enzyme - Tden) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        integration_den = odeint(PCR_reaction, values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

        concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den

        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tanneal * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tanneal * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tanneal * dS)

        dGs[3] = (Tm_enzyme - Tanneal) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        integration_anneal = odeint(PCR_reaction, integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)

        concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal

        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Text * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Text * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Text * dS)

        dGs[3] = (Tm_enzyme - Text) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)

        concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext

        values = integration_ext[-1]

    dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_cooling_down * dS)

    dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_cooling_down * dS)

    dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_cooling_down * dS)

    dGs[3] = (Tm_enzyme - T_cooling_down) * dS

    dGs = np.clip(dGs, a_min=None, a_max=1e+12)

    integration_cool = odeint(PCR_reaction, integration_ext[-1], time[(total * (i + 1) * steps) - 1: (total * (i + 1) * steps + t_cooling_down * steps)], args=(T_cooling_down, dGs), mxstep=5000000)

    concentration[(total * (i + 1) * steps) - 1: (total * (i + 1) * steps + t_cooling_down * steps)] = integration_cool[-1]

    values = integration_cool[-1]

    print("The concentration of the 17 species at the end of PCR integration is:", values)


    print("total conc", PCR_total_concentration(concentration, time))

    #after_nt_number = all_nucleotide(values)

    #all_umol_2 = u_molar_concentration(values)

    #print("The difference in nt number after", functions_plus_name[number], ":", before_nt_number - after_nt_number, "\n")

    #print("The difference in micromolar concentration after", functions_plus_name[number], ":", all_umol_1 - all_umol_2, "\n")

    #Plotting the concentrations over time

    plt.figure(1)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 18,  fontweight = 'bold', y = 0.95)

    plots = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(10):

        plt.subplot(2, 5, i+1)

        plt.gca().set_title(species[plots[i]], fontweight = 'bold')


        plt.plot(time, concentration[:, plots[i]])

        #plt.ylim([0, y_top_limit[i]])


        #plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time", fontsize = 12)
        plt.ylabel("Concentration", fontsize = 12)



    plt.subplots_adjust(wspace = 0.38)
    plt.show()

    return values



def PCR_total_concentration(all_concentration, time_vector):


    single_species_PCR = ["S1", "S2", "P1", "P2", "Q1", "Q2", "E"]


    concentration_PCR = np.zeros((all_concentration.shape[0], len(single_species_PCR)))         # the matrix s length is all time points

    indexes_of_species = [[] for i in range(len(single_species_PCR))]


    for x in range(len(single_species_PCR)):


        for i in species:

            if single_species_PCR[x] in i:


                indexes_of_species[x].append(species.index(i))


    print(indexes_of_species)




    for x in range(all_concentration.shape[0]):


        for i in range(len(indexes_of_species)):

            summary_concentration = 0


            for n in range(len(indexes_of_species[i])):


                summary_concentration = summary_concentration + all_concentration[x, indexes_of_species[i][n]]



            concentration_PCR[x, i] = summary_concentration




    plt.figure(1)

    plt.suptitle("Total concentrations of single species over time", FontSize= 16, FontWeight = "bold")


    #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(len(single_species_PCR)):

        plt.subplot(2, 4, i+1)


        plt.plot(time_vector, concentration_PCR[:, i], label = single_species_PCR[i])

        plt.gca().set_title(single_species_PCR[i])

        #plt.ylim([0, y_top_limit[i]])


        #plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time (s) ",  FontSize= 13, FontWeight = "bold")
        plt.ylabel("Total concentration (uM)",  FontSize= 13, FontWeight = "bold")
    #plt.tight_layout()




    plt.figure(2)

    plt.suptitle("Total concentrations of 4 single species over time\nwith accurate primer binding site" , fontsize=18,  fontweight= 'bold', y=0.95)

    highlighted_species = [0, 3, 5, 6]            #["S1", "S2", "P1", "P2", "Q1", "Q2", "E"]


    #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(len(highlighted_species)):

        plt.subplot(2, 2, i+1)

        plt.gca().set_title(single_species_PCR[highlighted_species[i]], fontweight = 'bold')


        plt.plot(time_vector, concentration_PCR[:, highlighted_species[i]])  # label = single_species_PCR[highlighted_species[i]]

        #plt.ylim([0, y_top_limit[i]])


        #plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time (s) ",  FontSize= 13)
        plt.ylabel("Total concentration (uM)",  FontSize= 13)
        #plt.legend(loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))




    plt.subplots_adjust(hspace = 0.3)

    plt.show()









    return concentration_PCR














if __name__ == '__main__':


    #print("dh", dH)
    # print(dH_check)
    print(Tm_S1S2)
    #print(Tm_enzyme2)




    #dS_change_4_species2(overall_concentration)


    #dS_change(values)


    #print(u_molar_concentration(values))

    #print(polymerase_nt_per_s(celsius_to_Kelvin(95)))

    #S1S2_dS(values)

    #dS1S2_dS_Tden(values)

    #print(taq_denaturation(values, Tanneal))


    #print(values)

    #individual_integration(values)

    #iterative_integration(values, 7)




    #print(len(new_species))

    PCR_integration(values)

    #print(cycle1)

    #print(values)

    #print(u_molar_concentration(values))
