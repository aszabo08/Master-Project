"""
This application visualizes primer binding and misbinding to a substrate, creating substrate-primer products.
"S" stands for the concentration of the substrate, "P" and "SP" for the concentration of the primer and the bound state of the substrate and primer.
The concentrations for all the above mentioned species are given in mole/liter.
The k parameters (k1, k2, k3) are rate constants in binding and misbinding processes.
"""


from scipy.integrate import odeint

import matplotlib.pyplot as plt

import numpy as np


#import extra


species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]


amplicon_length = 1000
primer_length = 15

dG_den = -2573.5784  # using RNAcofold with the two complementary sequences: dG = -615.10 kcal/mol which is equivalent with -2573.5784 kJ/mol




def denaturation(values, t, T):


    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]


    kf1 = 1
    kr1 = kf1 * np.exp(dG/(R*T))



    rate_den = -kf1 * S1S2 + kr1 * S1 * S2

    y = np.zeros(17)

    y[0] = rate_den
    y[1] = -rate_den
    y[2] = -rate_den


    return y



def primer_binding_1(values, t, T):


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

    kr2a = kf2 * np.exp(dG/(R*T))

    rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2

    kr2b = kf2 * np.exp(dG/(R*T))

    rate_S2P1_bind = kf2 * S2 * P1 - kr2b * S2P1


    y = np.zeros(17)


    y[1] = - rate_S1P2_bind
    y[2] = - rate_S2P1_bind
    y[3] = - rate_S2P1_bind
    y[4] = - rate_S1P2_bind
    y[5] = rate_S1P2_bind
    y[6] = rate_S2P1_bind


    return y



def primer_binding_2(values, t, T):


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

    kr5a = kf5 * np.exp(dG/(R*T))

    rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    kr5b = kf5 * np.exp(dG/(R*T))

    rate_S2Q1_bind = kf5 * S2 * Q1 - kr5b * S2Q1


    y = np.zeros(17)


    y[1] = - rate_S1Q2_bind
    y[2] = - rate_S2Q1_bind
    y[11] = - rate_S2Q1_bind
    y[12] = - rate_S1Q2_bind
    y[15] = rate_S1Q2_bind
    y[16] = rate_S2Q1_bind


    return y




def polymerase_binding_1(values, t, T):


    S1P2 = values[5]
    S2P1 = values[6]
    E = values[7]
    S1P2E = values[8]
    S2P1E = values[9]


    kf3 = 1


    kr3a = kf3 * np.exp(dG/(R*T))

    rate_poly_S1P2_bind = kf3 * S1P2 * E - kr3a * S1P2E


    kr3b = kf3 * np.exp(dG/(R*T))

    rate_poly_S2P1_bind = kf3 * S2P1 * E - kr3b * S2P1E


    y = np.zeros(17)

    y[5] = - rate_poly_S1P2_bind
    y[6] = - rate_poly_S2P1_bind
    y[7] = - rate_poly_S1P2_bind - rate_poly_S2P1_bind
    y[8] = rate_poly_S1P2_bind
    y[9] = rate_poly_S2P1_bind


    return y



def polymerase_binding_2(values, t, T):


    S1Q2 = values[15]
    S2Q1 = values[16]
    E = values[7]
    S1Q2E = values[13]
    S2Q1E = values[14]


    kf4 = 1


    kr4a = kf4 * np.exp(dG/(R*T))

    rate_poly_S1Q2_bind = kf4 * S1Q2 * E - kr4a * S1Q2E


    kr4b = kf4 * np.exp(dG/(R*T))

    rate_poly_S2Q1_bind = kf4 * S2Q1 * E - kr4b * S2Q1E


    y = np.zeros(17)

    y[15] = - rate_poly_S1Q2_bind
    y[16] = - rate_poly_S2Q1_bind
    y[7] = - rate_poly_S1Q2_bind - rate_poly_S2Q1_bind
    y[13] = rate_poly_S1Q2_bind
    y[14] = rate_poly_S2Q1_bind


    return y







def primer_ext_1(values, t, T):

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

    rate_ext_1 = ce / n * S1P2E * dNTP

    rate_ext_2 = ce / n * S2P1E * dNTP

    y = np.zeros(17)


    y[8] = - rate_ext_1    # concentration of S1P2E
    y[9] = - rate_ext_2   # concentration of S2P1E
    y[10] = - n * rate_ext_1 - n * rate_ext_2    # concentration of dNTP
    y[13] = rate_ext_1    # concentration of S1Q2E
    y[14] = rate_ext_2   # concentration of S2Q1E


    return y


def primer_ext_2(values, t, T):



    S1Q2E = values[13]
    S2Q1E = values[14]
    dNTP = values[10]

    ce_Q = 1000         # enzyme concentration

    # reaction: S1Q2E + extended_length * dNTP ---> S1S2 + E

    rate_ext_Q1 = ce_Q / extended_length * S1Q2E * dNTP


    # reaction: S2Q1E + extended_length * dNTP ---> S1S2 + E

    rate_ext_Q2 = ce_Q / extended_length * S2Q1E * dNTP


    y = np.zeros(17)

    y[0] = rate_ext_Q1 + rate_ext_Q2
    y[7] = rate_ext_Q1 + rate_ext_Q2
    y[10] = - n * rate_ext_Q1 - n * rate_ext_Q2    # concentration of dNTP
    y[13] = - rate_ext_Q1                           # concentration of S1Q2E
    y[14] = - rate_ext_Q2                           # concentration of S2Q1E




    return y


def PCR_reaction(values, t, T):

    return denaturation(values, t, T) + primer_binding_1(values, t, T) + primer_binding_2(values, t, T) + polymerase_binding_1(values, t, T) + polymerase_binding_2(values, t, T) + primer_ext_1(values, t, T) + primer_ext_2(values, t, T)




if __name__ == '__main__':

    values = [0 for i in range(17)]


    values[0] = float(input("Enter the concentration of plasmid (ng): "))


    values[3] = float(input("Enter the concentration of each primer (microL): "))


    values[4] = values[3]


    values[7] = float(input("Enter the concentration of polymerase (U): "))


    values[10] = float(input("Enter the concentration of each dNTP (microL): "))



    Tden_celsius, tden_string = input("Enter the temperature of denaturation (Celsius) and the length of it (second) ").split()

    # Converting the temperature from Celsius to Kelvin

    Tden = float(Tden_celsius) + 273.15

    tden = int(tden_string)


    Tanneal_celsius, tanneal_string = input("Enter the temperature of annealing (Celsius) and the length of it (second) ").split()

    # Converting the temperature from Celsius to Kelvin

    Tanneal = float(Tanneal_celsius) + 273.15

    tanneal = int(tanneal_string)



    Text_celsius, text_string = input("Enter the temperature of primer extension (Celsius) and the length of it (second) ").split()

    # Converting the temperature from Celsius to Kelvin

    Text = float(Text_celsius) + 273.15

    text = int(text_string)


    number_cycles = int(input("Enter the number of cycles "))

   # extra()


    total = tden + tanneal + text


    steps = 1


    time = np.linspace(0, number_cycles * (tden + tanneal + text), number_cycles * (tden + tanneal + text) * steps)  # for every second "steps" points are distinguished

    number_time_points = total * steps * number_cycles


    # 1 entropy unit = 4.184 J/ K m

    dS = - 52.7184           # entropy in J/ K m  ---->  we use dS = -12.6 entropy unit which is equal to −52.7184 J/ K m

    dH = dG_den + Tden * dS


    R = 8.314e-3        # Gas contant in  KJ / mol

    #kB = 1.38064852e-23  # Boltzmann constant"


    # Calculating the number of added nucleotides during the primer_ext_1

    # dG_ext = dH - Text * dS

    # 0.995 = np.exp(-dG_ext/(R*Text))

    dG_ext = - np.log(0.995) * R * Text

    # dG_ext = n * dH - Text * dS

    n = int(round((dG_ext + Text * dS) / dH))

    print("The number of added nucleotides to the primer during primer_ext_1:", n)


    extended_length = amplicon_length - (primer_length + n)




    concentration = np.empty((number_time_points, 17))


    for i in range(number_cycles):


        dG = dG_den

        integration_den = odeint(PCR_reaction, values, time[(total * i * steps) : ((total * i + tden) * steps)] , args=(Tden, ))

        concentration[(total * i * steps) : ((total * i + tden) * steps)] = integration_den


        dG = dH - Tanneal * dS

        integration_anneal = odeint(PCR_reaction, integration_den[-1], time[((total * i + tden) * steps)-1 : ((total * i + tden + tanneal) * steps)] , args=(Tanneal, ))

        concentration[((total * i + tden) * steps) - 1 : ((total * i + tden + tanneal) * steps)] = integration_anneal


        dG = dH - Text * dS

        integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) -1 : (total * (i + 1) * steps)], args=(Text, ))

        concentration[((total * i + tden + tanneal) * steps) -1 : (total * (i + 1) * steps)] = integration_ext


        values = integration_ext[-1]



    print("values", values)


    #Plotting the concentrartions over time in two figures

    plt.figure(1)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)

    for i in range(8):

        plt.subplot(2, 4, i+1)


        plt.plot(time, concentration[:, i])

        plt.legend([species[i]], loc='best', prop={'size':10})

        plt.xlabel("Time")
        plt.ylabel("Concentration")


    plt.figure(2)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)

    for i in range(9):


        plt.subplot(2, 5, i +1)


        plt.plot(time, concentration[:, i + 8])

        plt.legend([species[i + 8]], loc='best', prop={'size':10})

        plt.xlabel("Time")
        plt.ylabel("Concentration")

    plt.show()















    # when changing Tden to T in denaturation and other functions: TypeError: can't multiply sequence by non-int of type 'float'

    # to solve it I created a new T vector:

    # total = tden + tanneal + text
    #
    # Temp = np.empty(total * steps * number_cycles)
    #
    # for i in range(number_cycles):
    #
    #
    #
    #     Temp[(total * i * steps) : ((total * i + tden) * steps)] = Tden
    #
    #     Temp[((total * i + tden) * steps) : ((total * i + tden + tanneal) * steps)] = Tanneal
    #
    #     Temp[((total * i + tden + tanneal) * steps) : (total * (i + 1) * steps)] = Text

    # with this Temp vector I could multiply Temp with R in the backrate calculation

    # but then: ValueError: setting an array element with a sequence ( in the backrate calculations)

    # to solve it I changed the previously 1D array to 2D:  y = np.zeros((17, total * steps * number_cycles))

    # but then when running the code with odeint : RuntimeError: The array return by func must be one-dimensional, but got ndim=2.
