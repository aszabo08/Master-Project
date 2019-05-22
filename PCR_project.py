"""
This application visualizes primer binding and misbinding to a substrate, creating substrate-primer products.
"S" stands for the concentration of the substrate, "P" and "SP" for the concentration of the primer and the bound state of the substrate and primer.
The concentrations for all the above mentioned species are given in mole/liter.
The k parameters (k1, k2, k3) are rate constants in binding and misbinding processes.
"""


from scipy.integrate import odeint

#import matplotlib.pyplot as plt

import numpy as np






values = [0 for i in range(17)]

values[0] = 60      # concentration of plasmid (S1S2) in ng
values[1] = 0       # concentration of S1 in ng
values[2] = 0       # concentration of S2 in ng
values[3] = 9.2     # concentration of P1 in microL
values[4] = 9.2     # concentration of P2 in microL
values[5] = 0       # concentration of S1P2
values[6] = 0       # concentration of S2P1
values[7] = 8.2     # concentration of E
values[8] = 0       # concentration of S1P2E
values[9] = 0       # concentration of S2P1E
values[10] = 4.2    # concentration of dNTP
values[11] = 0      # concentration of Q1
values[12] = 0      # concentration of Q2
values[13] = 0      # concentration of S1Q2E
values[14] = 0      # concentration of S2Q1E
values[15] = 0      # concentration of S2Q1
values[16] = 0      # concentration of S1Q2

Tden = 363.15           # Kelvin which is equal to 90 degree Celsius
Tanneal = 345.15        # Kelvin which is equal to 72 degree Celsius
Text = 353.15           # Kelvin which is equal to 80 degree Celsius

tden = 2            # seconds
tanneal = 3         # seconds
text = 4            # seconds


total = tden + tanneal + text





number_cycles = 3

steps = 1

time = np.linspace(0, number_cycles * (tden + tanneal + text), number_cycles * (tden + tanneal + text) * steps)  # for every second "steps" points are distinguished


Temp = np.empty(total * steps * number_cycles)

number_time_points = total * steps * number_cycles


for i in range(number_cycles):



    Temp[(total * i * steps) : ((total * i + tden) * steps)] = Tden

    Temp[((total * i + tden) * steps) : ((total * i + tden + tanneal) * steps)] = Tanneal

    Temp[((total * i + tden + tanneal) * steps) : (total * (i + 1) * steps)] = Text



#T = np.array(number_cycles * ([([Tden] * tden * steps), ([Tanneal] * tanneal * steps), ([Text] * text * steps)]))

#Te = np.array([363.15, 363.15, 345.15, 345.15, 345.15, 353.15, 353.15, 353.15, 353.15])


#T = np.linspace()

dG_den = -2573.5784  # using RNAcofold with the two complementary sequences: dG = -615.10 kcal/mol which is equivalent with -2573.5784 kJ/mol

dG_bind1 = -1306.6632  # using RNAcofold with amplicon and primer: dG = -312.30 kcal/mol which is equivalent with -1306.6632 kJ/mol

dG_bind2 = -1294.948  # using RNAcofold with the complement of the amplicon and the complement of the primer: dG = -309.50 kcal/mol which is equivalent with -1294.948 kJ/mol

dG_bind_Qa = -800

dG_bind_Qb = -700

dGa = -500

dGb = -400

dG_pol_Qa = -1000

dG_pol_Qb = -1100


R = 8.314e-3        # Gas contant in  KJ / mol

#kB = 1.38064852e-23  # Boltzmann constant"



# 0.995 = np.exp(-(dG/kB*T))

# dG = - np.log(0.995) * kB * Text
#
# dG = n * dH - T * dS

n = 5  # number of bases in the initial extension



amplicon_length = 1000
primer_length = 15
extended_length = amplicon_length - (primer_length + n)




def denaturation(values, t, T):


    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]


    kf1 = 1
    kr1 = kf1 * np.exp(dG_den/(R*T))



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

    kr2a = kf2 * np.exp(dG_bind1/(R*T))

    rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2

    kr2b = kf2 * np.exp(dG_bind2/(R*T))

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

    kr5a = kf5 * np.exp(dG_bind_Qa/(R*T))

    rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    kr5b = kf5 * np.exp(dG_bind_Qb/(R*T))

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


    kr3a = kf3 * np.exp(dGa/(R*T))

    rate_poly_S1P2_bind = kf3 * S1P2 * E - kr3a * S1P2E


    kr3b = kf3 * np.exp(dGb/(R*T))

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


    kr4a = kf4 * np.exp(dG_pol_Qa/(R*T))

    rate_poly_S1Q2_bind = kf4 * S1Q2 * E - kr4a * S1Q2E


    kr4b = kf4 * np.exp(dG_pol_Qb/(R*T))

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


    print(Temp)


    concentration = np.empty((number_time_points, 17))




    for i in range(number_cycles):




        integration_den = odeint(PCR_reaction, values, time[(total * i * steps) : ((total * i + tden) * steps)] , args=(Tden, ))

        concentration[(total * i * steps) : ((total * i + tden) * steps)] = integration_den

        #concentration.append(integration_den)


        integration_anneal = odeint(PCR_reaction, integration_den[-1], time[(total * i + tden) * steps : ((total * i + tden + tanneal) * steps)] , args=(Tanneal, ))



        concentration[(total * i + tden) * steps : ((total * i + tden + tanneal) * steps)] = integration_anneal
        #concentration.append(integration_anneal)

        integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) : (total * (i + 1) * steps)], args=(Text, ))


        concentration[((total * i + tden + tanneal) * steps) : (total * (i + 1) * steps)] = integration_ext

        #concentration.append(integration_ext)

        values = integration_ext[-1]



    print("values", values)

    print(concentration)

    print(concentration[:, 0])

   # print("slice", concentration[0:][:, 0])

    #print("slice", concentration[1][: , 0])

    #print(len(concentration))










    # Plotting SP, SP_newk, SP1, SP1_newk in one figure


    # plt.figure(1)
    #
    # plt.suptitle("Change of concentration over time ", fontsize = 14)
    #
    #
    # # plot of SP
    #
    # plt.subplot(221)
    #
    # SP_plot = plt.plot(t, SP)
    #
    # plt.legend(["S", "P", "SP"], loc='upper left', prop={'size':10}, bbox_to_anchor=(1,1))
    #
    # plot_attributes(SP_plot)







    # "The species included in the reactions are the following: "

    # state vector:

    # values = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2, "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]

    #change = np.zeros((number_cycles, 17), dtype=np.int)

   #print(PCR_reaction(values, time, Te))



    #print(T)

    #print(R)



    # print(R*T)


    #
    #concentration = odeint(PCR_reaction, values, time, args=(Te, ))
    #
    #print(concentration)

    # for i in range(number_cycles):
    #
    #
    #     change[i] = odeint(PCR_reaction, time, T)
    #
    #     values = values + change[i]
    #
    # result = values




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
