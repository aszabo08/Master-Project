
"""
This application models the polymerase chain reaction with the opportunity of primer misbinding using the user' s input for initial concentrations (dsDNA, primer, dNTP, enzyme) and temperature settings. The amplified amount of dsDNA is displayed in micromolar and ng/ ul concentrations at the end of the process.

"""


from scipy.integrate import odeint

import matplotlib.pyplot as plt

import numpy as np


species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]


new_species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1", "S1M2", "S1M2E", "S1N2E", "S1L2", "L2", "S1N2", "L2P1", "L2P1E", "L2Q1E", "L1L2", "L2Q1", "L1", "L1P2", "L1P2E", "L1Q2E", "L1Q2"]


values = [0 for i in range(33)]


##########  Standard conditions  ##########

values[0] = 0.00015151515151515152        # 5 ng of 1000 bp dsDNA
values[3] = 0.5                           # concentration of P1 in uM
values[4] = 0.5                           # concentration of P2 in uM
values[7] = 0.2                           # concentration of E in uM
values[10] = 200                          # concentration of dNTP in uM

T_initial_den = 369.15                    # temperature of initial denaturation in K (96 °C)
Tden = 369.15                             # temperature of denaturation in K (96 °C)

Tanneal = 334.15                         # temperature of annealing in K (61 °C)
#Tanneal = 320.15
Text = 345.15                             # temperature of extension in K (72 °C)
T_cooling_down = 345.15                   # temperature of final extension in K (72 °C)

t_initial_den = 0                         # length of initial denaturation in seconds
tden = 10                                 # length of denaturation in seconds
tanneal = 10                              # length of annealing in seconds
text = 20                                 # length of extension in seconds
t_cooling_down =0                         # length of final extension in seconds

enzyme_type = 'taq'                       # either 'taq' or 'q5'

amplicon_length = 1000                    # bp
primer_length = 15                        # nt
n = 10                                    # nt

Tm_primer = 339.15                       # melting temperature of the primer in K (66 °C)
Tm_extended_primer = 346.15              # melting temperature of the extended primer in K (73 °C)

number_cycles = 32

#number_cycles = 3

##########  Experimental design for validation ##########

# values[0] = 0.00013730416992764068      # 10 ng of 2207 bp dsDNA
# values[3] = 0.25                        # concentration of P1 in uM
# values[4] = 0.25                        # concentration of P2 in uM
# values[7] = 0.16                        # concentration of E in uM
# values[10] = 200                        # concentration of dNTP in uM
#
# T_initial_den = 371.15                  # temperature of initial denaturation in K (98 °C)
# Tden = 371.15                           # temperature of denaturation in K (98 °C)
# # Tanneal = 345.15                      # temperature of annealing in K (72 °C)
# Tanneal = 338.15                        # temperature of annealing in K (65 °C)
# Text = 345.15                           # temperature of extension in K (72 °C)
# T_cooling_down = 345.15                 # temperature of final extension in K (72 °C)
#
# t_initial_den = 30                         # length of initial denaturation in seconds
# tden = 10                                  # length of denaturation in seconds
# # tanneal = 35                             # length of annealing in seconds
# # text = 35                                # length of extension in seconds
# tanneal = 20                               # length of annealing in seconds
# text = 50                                  # length of extension in seconds
# t_cooling_down = 300                       # length of final extension in seconds
#
# enzyme_type = 'q5'                         # either 'taq' or 'q5'
#
# amplicon_length = 2207                     # bp
# primer_length = 25                         # nt
# n = 10                                     # nt
#
# Tm_primer = 345.15                         # melting temperature of the primer in K (72 °C)
# Tm_extended_primer = 348.15                # melting temperature of the extended primer in K (75 °C)
#
# number_cycles = 32


##########  General settings ##########


functions_name = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer_ext_2"]

steps = 1


initial_fixed = values
total = tden + tanneal + text
extended_primer = primer_length + n
extended_length = amplicon_length - extended_primer
time = np.linspace(0, t_initial_den + number_cycles * (tden + tanneal + text) + t_cooling_down, number_cycles * (tden + tanneal + text) * steps + (t_cooling_down + t_initial_den) * steps)  # for every second "steps" points are distinguished
number_time_points = total * steps * number_cycles + (t_initial_den + t_cooling_down) * steps

dS = - 2                                    # entropy change, unit: kJ K^-1

#dS = - 4

#max_exponent = 15

max_exponent = 7

# min_clip = -1e+15
# max_clip = 1e+15
min_clip = -1e+12
max_clip = 1e+12

forward_rate = 1                                # uM/s
R = 8.314e-3                                    # Gas constant in kJ K^-1 mol^-1

Tmax = 373.15                                   # maximum melting temperature in K (100 °C)

K = (primer_length * extended_primer * (Tm_extended_primer - Tm_primer)) / (extended_primer * Tm_primer - primer_length * Tm_extended_primer)           # Michaelis constant

dH = (Tm_primer * dS * (primer_length + K)) / (Tmax * primer_length)                                 # enthalpy change, unit: kJ

Tm_S1S2 = (Tmax * amplicon_length * dH) / ((amplicon_length + K) * dS)

Tm_enzyme = 356.15                      # 85 °C



##########  Functions for conversions ##########

def uM_to_ng_per_ul(uM, bp) :


    """Converts micromolar concentration to ng/ul using 660 g/mol molecular mass for the DNA"""

    g_per_ul = uM * 1e-6 * 1e-6 * 660 * bp		# uM to M to mol/ul to g/ul

    return g_per_ul * 1e9


def celsius_to_Kelvin(x):

    """ Converts temperatures from degree Celisus to Kelvin"""

    return x + 273.15


##########  Rate limiting functions ##########

def clipping(kf, kr):

    """ By clipping function the ratio between the forward abd backward rate stays the same after the limitation"""

    return (kf / kr) * np.clip(kr, a_min= min_clip, a_max= max_clip), np.clip(kr, a_min= min_clip, a_max= max_clip)


def rate_clipping(x):

    """ The reaction rates are limited by a maximum and minimum value"""

    return np.clip(x, a_min= min_clip, a_max= max_clip)


def exponent_clipping(x):

    """ The maximum amount of an exponent is limited by this functions"""

    return np.clip(x, a_min= None, a_max= max_exponent)


##########  Additional functions ##########

def u_molar_concentration(values):

    """ This functions calculates the total micromolar concentration of the system"""

    length = [2 * amplicon_length, amplicon_length, amplicon_length, primer_length, primer_length, amplicon_length + primer_length, amplicon_length + primer_length, 0, amplicon_length + primer_length, amplicon_length + primer_length, 1, extended_primer, extended_primer, amplicon_length + extended_primer,  amplicon_length + extended_primer, amplicon_length + extended_primer, amplicon_length + extended_primer, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    return sum(values[i] * length[i] for i in range(len(values)))


def polymerase_nt_per_s(temperature):

    """ This function calculates the speed of DNA ploymerase, temperature is given in Kelvin"""


    temperature_scale = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(95), (95 - 0 + 1))  # x coordinates

    q5_temperature_data = [celsius_to_Kelvin(25), celsius_to_Kelvin(30), celsius_to_Kelvin(45), celsius_to_Kelvin(60), celsius_to_Kelvin(70), celsius_to_Kelvin(75), celsius_to_Kelvin(80), celsius_to_Kelvin(85), celsius_to_Kelvin(90), celsius_to_Kelvin(93) ]  # temperatures with known polymerase rate

    q5_rate_data = [0,  0.005, 0.010,  0.075,  0.11,  0.15,  0.12, 0.10, 0.020, 0]           # taq polymerase rate at given temperatures

    #used speed for validation:
    #q5_rate_data = [0.4 * q5_rate_data[i] for i in range(len(q5_rate_data))]

    q5_temp_interpolar = np.interp(temperature_scale, q5_temperature_data, q5_rate_data)        # calculating the taq polymerase rate between 0 and 95 degrees

    taq_temperature_data = [celsius_to_Kelvin(25), celsius_to_Kelvin(30), celsius_to_Kelvin(45), celsius_to_Kelvin(60), celsius_to_Kelvin(70), celsius_to_Kelvin(75), celsius_to_Kelvin(80), celsius_to_Kelvin(85), celsius_to_Kelvin(90), celsius_to_Kelvin(93) ]

    taq_rate_data = [0.75 * q5_rate_data[i] for i in range(len(q5_rate_data))]

    taq_temp_interpolar = np.interp(temperature_scale, taq_temperature_data, taq_rate_data)        # calculating the taq polymerase rate between 0 and 90 degrees

    # PLot of the enzyme's incorporation speed
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


        if ((temperature >= celsius_to_Kelvin(25)) and (temperature <= celsius_to_Kelvin(93))):

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




##########  Modelled interactions of the PCR ##########




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

    S1 = values[1]
    S2 = values[2]
    P1 = values[3]
    P2 = values[4]
    S1P2 = values[5]
    S2P1 = values[6]


    kf2 = forward_rate

    exponent_2a = exponent_clipping(dGs[1]/(R*T))

    kr2a = kf2 * np.exp(exponent_2a)

    kf2, kr2a = clipping(kf2, kr2a)

    rate_S1P2_bind = rate_clipping(kf2 * S1 * P2 - kr2a * S1P2)

    exponent_2b = exponent_clipping(dGs[1]/(R*T))

    kr2b = kf2 * np.exp(exponent_2b)

    kf2, kr2b = clipping(kf2, kr2b)

    rate_S2P1_bind = rate_clipping(kf2 * S2 * P1 - kr2b * S2P1)

    y = np.zeros(33)


    y[1] = - rate_S1P2_bind
    y[2] = - rate_S2P1_bind
    y[3] = - rate_S2P1_bind
    y[4] = - rate_S1P2_bind
    y[5] = rate_S1P2_bind
    y[6] = rate_S2P1_bind

    return y



def primer_binding_2(values, t, T, dGs):


    S1 = values[1]
    S2 = values[2]
    Q1 = values[11]
    Q2 = values[12]
    S1Q2 = values[15]
    S2Q1 = values[16]

    kf5 = forward_rate

    exponent_5a = exponent_clipping(dGs[2]/(R*T))

    kr5a = kf5 * np.exp(exponent_5a)

    kf5, kr5a = clipping(kf5, kr5a)

    rate_S1Q2_bind = rate_clipping(kf5 * S1 * Q2 - kr5a * S1Q2)

    exponent_5b = exponent_clipping(dGs[2]/(R*T))

    kr5b = kf5 * np.exp(exponent_5b)

    kf5, kr5b = clipping(kf5, kr5b)

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

    exponent_3a = exponent_clipping(dGs[3]/(R*T))

    kr3a = kf3 * np.exp(exponent_3a)

    kf3, kr3a = clipping(kf3, kr3a)
    rate_poly_S1P2_bind = rate_clipping(kf3 * S1P2 * E - kr3a * S1P2E)

    exponent_3b = exponent_clipping(dGs[3]/(R*T))

    kr3b = kf3 * np.exp(exponent_3b)

    kf3, kr3b = clipping(kf3, kr3b)

    rate_poly_S2P1_bind = rate_clipping(kf3 * S2P1 * E - kr3b * S2P1E)

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

    exponent_4a = exponent_clipping(dGs[3]/(R*T))

    kr4a = kf4 * np.exp(exponent_4a)

    kf4, kr4a = clipping(kf4, kr4a)

    rate_poly_S1Q2_bind = rate_clipping(kf4 * S1Q2 * E - kr4a * S1Q2E)

    exponent_4b = exponent_clipping(dGs[3]/(R*T))

    kr4b = kf4 * np.exp(exponent_4b)

    kf4, kr4b = clipping(kf4, kr4b)

    rate_poly_S2Q1_bind = rate_clipping(kf4 * S2Q1 * E - kr4b * S2Q1E)

    enzyme = rate_clipping(- rate_poly_S1Q2_bind - rate_poly_S2Q1_bind)


    y = np.zeros(33)

    y[15] = - rate_poly_S1Q2_bind
    y[16] = - rate_poly_S2Q1_bind
    y[7] = enzyme
    y[13] = rate_poly_S1Q2_bind
    y[14] = rate_poly_S2Q1_bind


    return y



def primer_ext_1(values, t, T, dGs):


    S1P2E = values[8]
    S2P1E = values[9]
    dNTP = values[10]


    #   S1P2E + n * dNTP -> S1Q2E
    #   S2P1E + n * dNTP -> S2Q1E

    ce = polymerase_nt_per_s(T)

    rate_ext_1 = rate_clipping((ce / n) * S1P2E * dNTP)

    rate_ext_2 = rate_clipping((ce / n) * S2P1E * dNTP)

    nucleotide = rate_clipping(- n * rate_ext_1 - n * rate_ext_2)

    y = np.zeros(33)

    y[8] = - rate_ext_1  # concentration of S1P2E
    y[9] = - rate_ext_2  # concentration of S2P1E
    y[10] = nucleotide  # concentration of dNTP
    y[13] = rate_ext_1  # concentration of S1Q2E
    y[14] = rate_ext_2  # concentration of S2Q1E

    return y



def primer_ext_2(values, t, T, dGs):


    S1Q2E = values[13]
    S2Q1E = values[14]
    dNTP = values[10]


    ce_Q = polymerase_nt_per_s(T)

    # reaction: S1Q2E + extended_length * dNTP ---> S1S2 + E

    rate_ext_Q1 = rate_clipping((ce_Q / extended_length) * S1Q2E * dNTP)

    # reaction: S2Q1E + extended_length * dNTP ---> S1S2 + E

    rate_ext_Q2 = rate_clipping((ce_Q / extended_length) * S2Q1E * dNTP)

    nucleotide_Q = rate_clipping(- extended_length * rate_ext_Q1 - extended_length * rate_ext_Q2)

    product = rate_clipping(rate_ext_Q1 + rate_ext_Q2)

    y = np.zeros(33)

    y[0] = product
    y[7] = product

    y[10] = nucleotide_Q  # concentration of dNTP
    y[13] = - rate_ext_Q1  # concentration of S1Q2E
    y[14] = - rate_ext_Q2  # concentration of S2Q1E

    return y



def enzyme_denaturation(values, t, T, dGs):

    E = values[7]

    y = np.zeros(33)

    rate = 0.0001           # uM ^-s

    taq_denaturation_T = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(110), 111)

    temperature_data = [celsius_to_Kelvin(0), celsius_to_Kelvin(90), celsius_to_Kelvin(110)]

    taq_denaturation = [0, 0, rate * E]

    # PLot of enzyme denaturation

    temp_interpolar = np.interp(taq_denaturation_T, temperature_data, taq_denaturation)        # calculating the taq polymerase rate between 0 and 90 degrees

    # plt.title("The level of denaturation of DNA polymerase at different temperatures")
    #
    # plt.plot(taq_denaturation_T, temp_interpolar)
    #
    # plt.xlabel("Temperature (Kelvin)")
    # plt.ylabel("Level of DNA polymerase denaturation (uM) ")
    #
    # plt.show()

    if T > celsius_to_Kelvin(90):

        y[7] = - rate * E

    return y




def PCR_reaction(values, t, T, dGs):

    """ The summary of every interaction in the PCR"""

    summary = denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + primer_binding_2(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + primer_ext_2(values, t, T, dGs) + enzyme_denaturation(values, t, T, dGs)

    summary = np.clip(summary, a_min= min_clip, a_max= max_clip)

    return summary




def enzyme_dissociation(values, t, T, dGs):

    """ The summary of two functions: primer binding 1 and polymerase binding 1"""

    return primer_binding_1(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs)



##########  Assessing the value of dS ##########



low_concentration = [0 for i in range(33)]

low_concentration[0] = 0.00015151515151515152                       # 5 ng of 1000bp long DNA

low_concentration[5], low_concentration[6] = 0.0001, 0.0001

low_concentration[8], low_concentration[9] = 0.0001, 0.0001

low_concentration[15], low_concentration[16] = 0.0001,  0.0001


high_concentration = [0 for i in range(33)]

high_concentration[0] = 0.0175

high_concentration[5], high_concentration[6] = 0.0175, 0.0175

high_concentration[8], high_concentration[9] = 0.007, 0.007

high_concentration[15], high_concentration[16] = 0.014, 0.014

overall_concentration = [low_concentration, high_concentration]

initial_overall = overall_concentration






def dS_change_4_species2(overall_concentration):

    """ This function shows the dehybridisation behaviour of the chosen complexes around their melting temperatures"""


    functions = [denaturation, primer_binding_1, primer_binding_2, enzyme_dissociation]

    dGs = [0 for x in range(4)]

    var_dS = [-10, -5, -1, -0.5]

    time_all = np.linspace(0, tden + tanneal + text, tden + tanneal + text)


    temperature_scale = np.linspace(celsius_to_Kelvin(0), celsius_to_Kelvin(110), 111)     # Temperature scale between 0 and 110 degree Celsius

    var_dS_S1S2_result = np.zeros((len(overall_concentration), len(var_dS), len(temperature_scale)))

    var_dS_primer_result = np.zeros((len(overall_concentration), len(var_dS), len(temperature_scale)))

    var_dS_ext_primer_result = np.zeros((len(overall_concentration), len(var_dS), len(temperature_scale)))

    var_dS_enzyme_result = np.zeros((len(overall_concentration), len(var_dS), len(temperature_scale)))


    for e in range(len(overall_concentration)):


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

    # Plots of the dissociation curves

    size = 20

    l_size = 18

    plt.figure(1)

    plt.title("Low (-) and high (:) S1S2 concentration with different dS values", FontSize= 22, FontWeight = "bold", position=(0.5, 1.05))

    style_curve = ['-',':']

    colour_curve = ['C1', 'C2', 'C3','C4']


    for i in range(len(overall_concentration)):


        for x in range(len(var_dS)):


            plt.plot(temperature_scale, var_dS_S1S2_result[i, x, ]/initial_overall[i][0]*100, style_curve[i], color = colour_curve[x], label = "dS = " + str(var_dS[x]))     #/initial_overall[e][0]*100



    plt.axvline(x=Tm_S1S2, color = 'black', linestyle= "--", label = "Tm_S1S2")


    plt.xlim(310, 360)

    plt.xlabel("Temperature (K)", FontSize= size, FontWeight = "bold")

    plt.tick_params(labelsize = 16)

    plt.ylabel("Percentage of S1S2 concentration \n at the end of denaturation \n compared to the initial S1S2 concentration", FontSize= size, FontWeight = "bold")

    plt.legend(loc='upper left', prop={'size':l_size}, bbox_to_anchor=(1,1))

    plt.figure(2)

    plt.title("Low (-) and high (:) S1P2 concentration with different dS values", FontSize= 22, FontWeight = "bold", position=(0.5, 1.05))


    for i in range(len(overall_concentration)):


        for x in range(len(var_dS)):


            plt.plot(temperature_scale, var_dS_primer_result[i, x, ]/initial_overall[i][5]*100, style_curve[i], color = colour_curve[x], label = "dS = " + str(var_dS[x]))          #/initial_overall[e][5]*100



    plt.axvline(x=Tm_primer, color = 'black', linestyle= "--", label = "Tm_primer")

    plt.xlim(300,350)

    plt.xlabel("Temperature (K)", FontSize= size, FontWeight = "bold")

    plt.ylabel("Percentage of S1P2 concentration \n at the end of primer binding 1 reaction\n compared to the initial S1P2 concentration", FontSize= size, FontWeight = "bold")

    plt.legend(loc='upper left', prop={'size':l_size}, bbox_to_anchor=(1,1))


    plt.tick_params(labelsize = 16)


    plt.figure(3)

    plt.title("Low (-) and high (:) S1Q2 concentration with different dS values", FontSize= 22, FontWeight = "bold", position=(0.5, 1.05))


    for i in range(len(overall_concentration)):


        for x in range(len(var_dS)):

            plt.plot(temperature_scale, var_dS_ext_primer_result[i, x, ]/initial_overall[i][15]*100, style_curve[i], color = colour_curve[x], label = "dS = " + str(var_dS[x]))      # /initial_overall[e][15]*100


    plt.axvline(x=Tm_extended_primer, color = 'black', linestyle= "--", label = "Tm_ext_p")

    plt.xlim(300,360)

    plt.xlabel("Temperature (K)", FontSize= size, FontWeight = "bold")

    plt.ylabel("Percentage of S1Q2 concentration \n at the end of primer binding 2 reaction\n compared to the initial S1Q2 concentration", FontSize= size, FontWeight = "bold")

    plt.legend(loc='upper left', prop={'size':l_size}, bbox_to_anchor=(1,1))

    plt.tick_params(labelsize = 16)

    plt.figure(4)

    plt.title("Low (-) and high (:) S1P2E concentration with different dS values", FontSize= 22, FontWeight = "bold", position=(0.5, 1.05))


    for i in range(len(overall_concentration)):


        for x in range(len(var_dS)):


            plt.plot(temperature_scale, var_dS_enzyme_result[i, x, ]/initial_overall[i][8]*100, style_curve[i], color = colour_curve[x], label = "dS = " + str(var_dS[x]))           # /initial_overall[e][8]*100



    plt.axvline(x=Tm_enzyme, color = 'black', linestyle= "--", label = "Tm_enzyme")

    plt.xlim(300,360)

    plt.xlabel("Temperature (K)", FontSize= size, FontWeight = "bold")


    plt.ylabel("Percentage of S1P2E concentration \n at the end of enzyme_dissociation reaction\n compared to the initial S1P2E concentration", FontSize= size, FontWeight = "bold")


    plt.legend(loc='upper left', prop={'size':l_size}, bbox_to_anchor=(1,1))

    plt.tick_params(labelsize = 16)


    plt.show()


##########  Integration functions ##########



def individual_integration(values):

    """ The interactions of the PCR are run individially, one after the other

    The resulted plot informs about the total micromolar concentration difference"""

    functions = [denaturation, primer_binding_1, polymerase_binding_1, primer_ext_1, polymerase_binding_2, primer_binding_2, primer_ext_2]

    functions_name = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer ext_2"]


    concentration = np.empty((number_time_points, 33))

    dGs = [0 for x in range(4)]

    difference_umol = []

    for n in range(len(functions)):


        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_initial_den * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_initial_den * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_initial_den * dS)

        dGs[3] = (Tm_enzyme - T_initial_den) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        # when initial denaturation and final elongation are used, uncomment the below sections


        # integration_initial_den = odeint(PCR_reaction, values, time[0: t_initial_den * steps], args=(T_initial_den, dGs), mxstep=5000000)
        #
        # concentration[0: t_initial_den * steps] = integration_initial_den
        #
        # values = integration_initial_den[-1]

        all_umol_1 = u_molar_concentration(values)

        for i in range(number_cycles):

            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tden * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tden * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tden * dS)

            dGs[3] = (Tm_enzyme - Tden) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_den = odeint(PCR_reaction, values, time[(t_initial_den * steps -1 + total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)
            #
            # concentration[t_initial_den * steps -1 + (total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)] = integration_den


            integration_den = odeint(functions[n], values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

            concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den


            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tanneal * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tanneal * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tanneal * dS)

            dGs[3] = (Tm_enzyme - Tanneal) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_anneal = odeint(PCR_reaction, integration_den[-1], time[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)
            #
            # concentration[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)] = integration_anneal


            integration_anneal = odeint(functions[n], integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)

            concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal


            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Text * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Text * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Text * dS)

            dGs[3] = (Tm_enzyme - Text) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)
            #
            # concentration[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1:t_initial_den * steps +  (total * (i + 1) * steps)] = integration_ext


            integration_ext = odeint(functions[n], integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)

            concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


            values = integration_ext[-1]

        # dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_cooling_down * dS)
        #
        # dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_cooling_down * dS)
        #
        # dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_cooling_down * dS)
        #
        # dGs[3] = (Tm_enzyme - T_cooling_down) * dS
        #
        # dGs = np.clip(dGs, a_min=None, a_max=1e+12)
        #
        # integration_cool = odeint(PCR_reaction, integration_ext[-1], time[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)], args=(T_cooling_down, dGs), mxstep=5000000)
        #
        # concentration[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)] = integration_cool[-1]
        #
        # values = integration_cool[-1]


        all_umol_2 = u_molar_concentration(values)

        difference_umol.append(all_umol_1 - all_umol_2)


        print("The concentration of the 17 species after", functions_name[n], "is added to the process individually:", values)



        #print("The difference in micromol after", functions_name[n], ":", all_umol_1 - all_umol_2, "\n")

        #
        # plt.figure(1)
        #
        # plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)
        #
        # plots = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]
        #
        #
        # for i in range(10):
        #
        #     plt.subplot(2, 5, i+1)
        #
        #
        #     plt.plot(time, concentration[:, plots[i]])
        #
        #
        #     plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})
        #
        #     plt.xlabel("Time")
        #     plt.ylabel("Concentration")



        #plt.show()

    print("diff", difference_umol)

    plt.figure(2)

    plt.title("Difference of the overall micromolar concentration after each function is added individually", FontSize= 22, FontWeight = "bold", position=(0.5, 1.05))

    plt.plot(range(len(functions)), difference_umol)

    plt.tick_params(labelsize = 18)

    plt.xticks(range(len(functions)), functions_name)

    plt.xlabel("Name of the individually added function", FontSize= 20, FontWeight = "bold")
    plt.ylabel("Difference in micromolar concentration", FontSize= 20, FontWeight = "bold")

    plt.show()


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

    """ The interactions of the PCR are run iteratively,
    the maximum number of interaction is set by "number"
    (for the whole PCR 7 is the maximum number)

    The resulted plot informs about the total micromolar concentration difference"""



    all_summaries = [a,b,c, d,e,f,g]

    difference_umol = []

    for n in range(number):


        concentration = np.empty((number_time_points, 33))

        dGs = [0 for x in range(4)]

        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_initial_den * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_initial_den * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_initial_den * dS)

        dGs[3] = (Tm_enzyme - T_initial_den) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)


        # integration_initial_den = odeint(PCR_reaction, values, time[0: t_initial_den * steps], args=(T_initial_den, dGs), mxstep=5000000)
        #
        # concentration[0: t_initial_den * steps] = integration_initial_den
        #
        # values = integration_initial_den[-1]

        all_umol_1 = u_molar_concentration(values)

        for i in range(number_cycles):

            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tden * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tden * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tden * dS)

            dGs[3] = (Tm_enzyme - Tden) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_den = odeint(PCR_reaction, values, time[(t_initial_den * steps -1 + total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)
            #
            # concentration[t_initial_den * steps -1 + (total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)] = integration_den


            integration_den = odeint(all_summaries[n], values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

            concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den


            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tanneal * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tanneal * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tanneal * dS)

            dGs[3] = (Tm_enzyme - Tanneal) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_anneal = odeint(PCR_reaction, integration_den[-1], time[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)
            #
            # concentration[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)] = integration_anneal


            integration_anneal = odeint(all_summaries[n], integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)

            concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal



            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Text * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Text * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Text * dS)

            dGs[3] = (Tm_enzyme - Text) * dS

            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)
            #
            # concentration[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1:t_initial_den * steps +  (total * (i + 1) * steps)] = integration_ext


            integration_ext = odeint(all_summaries[n], integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)

            concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


            values = integration_ext[-1]

        # dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_cooling_down * dS)
        #
        # dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_cooling_down * dS)
        #
        # dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_cooling_down * dS)
        #
        # dGs[3] = (Tm_enzyme - T_cooling_down) * dS
        #
        # dGs = np.clip(dGs, a_min=None, a_max=1e+12)
        #
        # integration_cool = odeint(PCR_reaction, integration_ext[-1], time[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)], args=(T_cooling_down, dGs), mxstep=5000000)
        #
        # concentration[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)] = integration_cool[-1]
        #
        # values = integration_cool[-1]


        all_umol_2 = u_molar_concentration(values)

        difference_umol.append(all_umol_1 - all_umol_2)

        print("The concentration of the 17 species after", functions_name[n],  "is added to the process iteratively:", values)


        #print("The difference in micromol after", functions_name[n], ":", all_umol_1 - all_umol_2, "\n")


        # plt.figure(1)
        #
        # plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)
        #
        # plots = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]
        #
        #
        # for i in range(10):
        #
        #     plt.subplot(2, 5, i+1)
        #
        #
        #     plt.plot(time, concentration[:, plots[i]])
        #
        #
        #     plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})
        #
        #     plt.xlabel("Time")
        #     plt.ylabel("Concentration")
        #
        #
        # plt.show()
        #
        # values = initial_fixed


    plt.figure(2)

    plt.title("Difference of micromolar concentration after each function is added iteratively",  FontSize= 18, FontWeight = "bold", position=(0.5, 1.05))

    plt.plot(range(number), difference_umol)

    plt.tick_params(labelsize = 14)

    plt.xticks(range(number), functions_name)

    plt.xlabel("Name of the iteratively added function", FontSize= 14, FontWeight = "bold")
    plt.ylabel("Difference in micromolar concentration", FontSize= 14, FontWeight = "bold")

    plt.show()


    return integration_ext[-1]



def PCR_integration(values):

    """ Integrating the PCR reaction over t time"""

    concentration = np.empty((number_time_points, 33))          # len new species!

    dGs = [0 for x in range(4)]

    #all_umol_1 = u_molar_concentration(values)

    dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_initial_den * dS)

    dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_initial_den * dS)

    dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_initial_den * dS)

    dGs[3] = (Tm_enzyme - T_initial_den) * dS

    dGs = np.clip(dGs, a_min=None, a_max=1e+12)



    # integration_initial_den = odeint(PCR_reaction, values, time[0: t_initial_den * steps], args=(T_initial_den, dGs), mxstep=5000000)
    #
    # concentration[0: t_initial_den * steps] = integration_initial_den
    #
    # values = integration_initial_den[-1]


    for i in range(number_cycles):

        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tden * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tden * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tden * dS)

        dGs[3] = (Tm_enzyme - Tden) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        # integration_den = odeint(PCR_reaction, values, time[(t_initial_den * steps -1 + total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)
        #
        # concentration[t_initial_den * steps -1 + (total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)] = integration_den


        integration_den = odeint(PCR_reaction, values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

        concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den


        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tanneal * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tanneal * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tanneal * dS)

        dGs[3] = (Tm_enzyme - Tanneal) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        # integration_anneal = odeint(PCR_reaction, integration_den[-1], time[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)
        #
        # concentration[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)] = integration_anneal


        integration_anneal = odeint(PCR_reaction, integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)

        concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal

        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Text * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Text * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Text * dS)

        dGs[3] = (Tm_enzyme - Text) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        # integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)
        #
        # concentration[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1:t_initial_den * steps +  (total * (i + 1) * steps)] = integration_ext


        integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)

        concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext

        values = integration_ext[-1]

    # dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_cooling_down * dS)
    #
    # dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_cooling_down * dS)
    #
    # dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_cooling_down * dS)
    #
    # dGs[3] = (Tm_enzyme - T_cooling_down) * dS
    #
    # dGs = np.clip(dGs, a_min=None, a_max=1e+12)
    #
    # integration_cool = odeint(PCR_reaction, integration_ext[-1], time[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)], args=(T_cooling_down, dGs), mxstep=5000000)
    #
    # concentration[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)] = integration_cool[-1]
    #
    # values = integration_cool[-1]

    print("The concentration of the 17 species at the end of PCR integration is:", values)

    print("The concentration of S1S2 in ng/ul:", uM_to_ng_per_ul(values[0], amplicon_length))


    #print("total conc", PCR_total_concentration(concentration, time))

    #all_umol_2 = u_molar_concentration(values)


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

        plt.xlabel("Time (s)", fontsize = 12)
        plt.ylabel("Concentration (uM)", fontsize = 12)

    plt.subplots_adjust(wspace = 0.38)
    plt.show()

    return values



def PCR_total_concentration(all_concentration, time_vector):

    """ This function calculates the total concentrations of the main species"""


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

        plt.xlabel("Time (s) ",  FontSize= 13, FontWeight = "bold")
        plt.ylabel("Total concentration (uM)",  FontSize= 13, FontWeight = "bold")



    plt.figure(2)

    plt.suptitle("The change of total concentrations of 4 single species\nwith accurate primer binding site" , fontsize=22,  fontweight= 'bold')

    highlighted_species = [0, 3, 5, 6]            #["S1", "S2", "P1", "P2", "Q1", "Q2", "E"]


    #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(len(highlighted_species)):

        plt.subplot(2, 2, i+1)

        plt.tick_params(labelsize = 16)

        plt.gca().set_title(single_species_PCR[highlighted_species[i]], fontweight = 'bold', fontsize = 18)


        plt.plot(time_vector, concentration_PCR[:, highlighted_species[i]])  # label = single_species_PCR[highlighted_species[i]]

        #plt.ylim([0, y_top_limit[i]])


        plt.xlabel("Time (s) ",  FontSize= 18)
        plt.ylabel("Total concentration (uM)",  FontSize= 18)

    plt.subplots_adjust(hspace = 0.3)

    plt.show()



    return concentration_PCR




if __name__ == '__main__':












    PCR_integration(values)


