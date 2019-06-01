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

values[0] = 1.5385e-15      # concentration of plasmid (S1S2) in mol/ul
values[1] = 0       # concentration of S1 in ng
values[2] = 0       # concentration of S2 in ng
values[3] = 1.5385e-11    # concentration of P1 in mol/ul
values[4] = 1.5385e-11   # concentration of P2 in mol/ul
values[5] = 0       # concentration of S1P2
values[6] = 0       # concentration of S2P1
values[7] = 1.6e-4  # concentration of E in mol/microl
values[8] = 0       # concentration of S1P2E
values[9] = 0       # concentration of S2P1E
values[10] = 1.5385e-11  # concentration of dNTP in mol/ul
values[11] = 0      # concentration of Q1
values[12] = 0      # concentration of Q2
values[13] = 0      # concentration of S1Q2E
values[14] = 0      # concentration of S2Q1E
values[15] = 0      # concentration of S2Q1
values[16] = 0      # concentration of S1Q2

#Tden = 363.15           # Kelvin which is equal to 90 degree Celsius

Tden = 340

#Tanneal = 345.15        # Kelvin which is equal to 72 degree Celsius

Tanneal = 300


#Text = 353.15           # Kelvin which is equal to 80 degree Celsius

Text = 350

tden = 2            # seconds
tanneal = 4      # seconds
text = 4            # seconds


total = tden + tanneal + text

number_cycles = 2

steps = 1


amplicon_length = 1000
primer_length = 15

n = 10

extended_primer = primer_length + n

extended_length = amplicon_length - extended_primer

time = np.linspace(0, number_cycles * (tden + tanneal + text), number_cycles * (tden + tanneal + text) * steps)  # for every second "steps" points are distinguished

number_time_points = total * steps * number_cycles

# dGs = [dG_S1S2, dG_S1P2_S2P1, dG_S1Q2_S2Q1, dG_enzyme]


#T = np.array(number_cycles * ([([Tden] * tden * steps), ([Tanneal] * tanneal * steps), ([Text] * text * steps)]))


# 1 entropy unit = 4.184 J/ K m

dS = - 52.7184          # entropy in J/ K m  ---->  we use dS = -12.6 entropy unit which is equal to −52.7184 J/ K m which is equal to - 0.0527184 kJ / K m

#dG_den = -257.5784  # using RNAcofold with the two complementary sequences: dG = -615.10 kcal/mol which is equivalent with -2573.5784 kJ/mol



# At melting temperature dG = 0

# dG = n * dH - Tm * dS    -----> n * dH = Tm * dS    -----> dH = 1/n * Tm * dS

Tm_S1S2 = 358.15        # in Kelvin which is equal to 85 degree Celsius

Tm_S1P2_S2P1 = 303.15        # in Kelvin which is equal to 30 degree Celsius

                            # Tm_S1P2 = Tm_S2P1



# At melting temperature dG = 0

# dH/n and dS are the same for the substrates, primer and extended primer

# amplicon_length * dH - Tm_S1S2 * dS = primer_length * dH - Tm_S1P2_S2P1 * dS
#
# amplicon_length * dH - primer_length * dH = Tm_S1S2 * dS - Tm_S1P2_S2P1 * dS
#
# (amplicon_length - primer_length) * dH = (Tm_S1S2 - Tm_S1P2_S2P1) * dS
#
# dH = ((Tm_S1S2 - Tm_S1P2_S2P1) * dS) / (amplicon_length - primer_length)
#
# dS = ((amplicon_length - primer_length) * dH) / (Tm_S1S2 - Tm_S1P2_S2P1)
#
# dS = (amplicon_length * dH) / Tm_S1S2 = (primer_length * dH) / Tm_S1P2_S2P1
#
#
#
# dH = Tm_S1S2 * dS / amplicon_length









Tm_S1Q2_S2Q1 = 323.15       # in Kelvin which is equal to 50 degree Celsius





dH = - 2.9437

dH_enzyme2 = 5 * dH


#
# dH_S1S2 = (Tm_S1S2 * dS) / amplicon_length              # J/ m
#
# dH_S1P2_S2P1 = (Tm_S1P2_S2P1 * dS) / primer_length
#
#
# dH_S1Q2_S2Q1 = (Tm_S1Q2_S2Q1 * dS) / extended_primer
#
# dH_enzyme = 30 * dH_S1P2_S2P1



# dG_S1S2 = amplicon_length * dH_S1S2 - T * dS          # J/ m
#
# dG_S1P2_S2P1 = primer_length * dH_S1P2_S2P1 - T * dS
#
# dG_S1Q2_S2Q1 = extended_primer * dH_S1Q2_S2Q1 - T * dS
#
# dG_enzyme = dH_enzyme - T * dS





# or Calculating the number of added nucleotides during the primer_ext_1
#
# # dG_ext = dH - Text * dS
#
# # 0.995 = np.exp(-dG_ext/(R*Text))
#
# dG_ext = - np.log(0.995) * R * T
#
#  dG_ext = n * dH - T * dS
#
# n = int(round((dG_ext + T * dS) / dH))









R = 8.314e-3        # Gas contant in  J / K mol

#kB = 1.38064852e-23  # Boltzmann constant"


#

print("The number of added nucleotides to the primer during primer_ext_1:", n)







def denaturation(values, t, T, dGs):


    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]

    print("S1S2, S1, S2", S1S2, S1, S2 )


    kf1 = 1


    exponent_1 = dGs[0]/(R*T)

    exponent_1 = np.clip(exponent_1, a_min= None,  a_max= 709.7827128933827 )

    kr1 = kf1 * np.exp(exponent_1)

    print("kr1 origin", kr1)

    kr1 = np.clip(kr1, a_min= -1e+14, a_max= 1e+14)

    print("kr1 second", kr1)


    rate_den = -kf1 * S1S2 + kr1 * S1 * S2

    rate_den = np.clip(rate_den,  a_min= -1e+14, a_max= 1e+14 )

    print("rateden", rate_den)

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

    print("exponent_2a origin", exponent_2a)

    exponent_2a = np.clip(exponent_2a, a_min= None, a_max= 709.7827128933827)

    print("exponent_2a second", exponent_2a)

    kr2a = kf2 * np.exp(exponent_2a)

    kr2a = np.clip(kr2a,  a_min= -1e+14, a_max= 1e+14 )



    print("kr2a second", kr2a)



    rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2


    rate_S1P2_bind = np.clip(rate_S1P2_bind,  a_min= -1e+14, a_max= 1e+14 )

    print("rate S1P1 bind", rate_S1P2_bind)

    exponent_2b = dGs[1]/(R*T)

    print("exponent_2b origin", exponent_2b)

    exponent_2b = np.clip(exponent_2b, a_min= None, a_max= 709.7827128933827)

    print("exponent_2b second", exponent_2b)

    kr2b = kf2 * np.exp(exponent_2b)


    kr2b = np.clip(kr2b,  a_min= -1e+14, a_max= 1e+14 )

    print("kr2b", kr2b)

    rate_S2P1_bind = kf2 * S2 * P1 - kr2b * S2P1


    rate_S2P1_bind = np.clip(rate_S2P1_bind,  a_min= -1e+14, a_max= 1e+14 )

    print("rate_s2P1", rate_S2P1_bind)


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

    exponent_5a = np.clip(exponent_5a, a_min= None, a_max= 709.7827128933827)

    print("exponent_5a", exponent_5a)

    kr5a = kf5 * np.exp(exponent_5a)

    kr5a = np.clip(kr5a,  a_min= -1e+14, a_max= 1e+14 )

    print("kr5a", kr5a)

    rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    rate_S1Q2_bind = np.clip(rate_S1Q2_bind,  a_min= -1e+14, a_max= 1e+14 )

    print("rate_S1Q2_bind", rate_S1Q2_bind)

    exponent_5b = dGs[2]/(R*T)

    exponent_5b = np.clip(exponent_5b, a_min= None, a_max= 709.7827128933827)


    print("exponent_5b", exponent_5b)

    kr5b = kf5 * np.exp(exponent_5b)


    kr5b = np.clip(kr5b,  a_min= -1e+14, a_max= 1e+14 )

    print("kr5b", kr5b)

    rate_S2Q1_bind = kf5 * S2 * Q1 - kr5b * S2Q1

    rate_S2Q1_bind = np.clip(rate_S2Q1_bind,  a_min= -1e+14, a_max= 1e+14 )

    print("rate_S2Q1_bind", rate_S2Q1_bind)

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

    exponent_3a = np.clip(exponent_3a, a_min= None, a_max= 709.7827128933827)

    kr3a = kf3 * np.exp(exponent_3a)

    kr3a = np.clip(kr3a,  a_min= -1e+14, a_max= 1e+14 )

    rate_poly_S1P2_bind = kf3 * S1P2 * E - kr3a * S1P2E

    rate_poly_S1P2_bind = np.clip(rate_poly_S1P2_bind,  a_min= -1e+14, a_max= 1e+14 )

    exponent_3b = dGs[3]/(R*T)

    exponent_3b = np.clip(exponent_3b, a_min= None, a_max= 709.7827128933827)


    kr3b = kf3 * np.exp(exponent_3b)

    kr3b = np.clip(kr3b,  a_min= -1e+14, a_max= 1e+14 )

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

    exponent_4a = np.clip(exponent_4a, a_min= None, a_max= 709.7827128933827)


    kr4a = kf4 * np.exp(exponent_4a)

    kr4a = np.clip(kr4a,  a_min= -1e+14, a_max= 1e+14)

    rate_poly_S1Q2_bind = kf4 * S1Q2 * E - kr4a * S1Q2E

    rate_poly_S1Q2_bind = np.clip(rate_poly_S1Q2_bind,  a_min= -1e+14, a_max= 1e+14)



    exponent_4b = dGs[3]/(R*T)

    exponent_4b = np.clip(exponent_4b, a_min= None, a_max= 709.7827128933827)

    print("exponent_4b", exponent_4b)


    kr4b = kf4 * np.exp(exponent_4b)

    kr4b = np.clip(kr4b,  a_min= -1e+14, a_max= 1e+14 )

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

    print("rate_exp_1", rate_ext_1)

    rate_ext_2 = ce / n * S2P1E * dNTP

    rate_ext_2 = np.clip(rate_ext_2,  a_min= -1e+14, a_max= 1e+14)

    print("rate_ext_2", rate_ext_2)

    nucleotide = - n * rate_ext_1 - n * rate_ext_2

    nucleotide = np.clip(nucleotide,  a_min= -1e+14, a_max= 1e+14)

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

    nucleotide_Q = - n * rate_ext_Q1 - n * rate_ext_Q2

    nucleotide_Q = np.clip(nucleotide_Q,  a_min= -1e+14, a_max= 1e+14)

    y = np.zeros(17)

    y[0] = rate_ext_Q1 + rate_ext_Q2
    y[7] = rate_ext_Q1 + rate_ext_Q2
    y[10] = nucleotide_Q   # concentration of dNTP
    y[13] = - rate_ext_Q1                           # concentration of S1Q2E
    y[14] = - rate_ext_Q2                           # concentration of S2Q1E




    return y


def PCR_reaction(values, t, T, dGs):

    summary = denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + primer_binding_2(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + primer_ext_2(values, t, T, dGs)

    summary = np.clip(summary,  a_min= -1e+14, a_max= 1e+14)

    return summary


if __name__ == '__main__':



    concentration = np.empty((number_time_points, 17))

    dGs = [0,0,0,0]




    for i in range(number_cycles):


        # dGs[0] = amplicon_length * dH_S1S2 - Tden * dS
        #
        # dGs[1] = primer_length * dH_S1P2_S2P1 - Tden * dS
        #
        # dGs[2] = extended_primer * dH_S1Q2_S2Q1 - Tden * dS
        #
        # dGs[3] = dH_enzyme - Tden * dS



        dGs[0] = amplicon_length * dH - Tden * dS

        dGs[1] = primer_length * dH - Tden * dS

        dGs[2] = extended_primer * dH - Tden * dS

        dGs[3] = dH_enzyme2 - Tden * dS

        dGs = np.clip(dGs, a_min= None, a_max= 1e+14 )

        print("dGs_den", dGs)



        integration_den = odeint(PCR_reaction, values, time[(total * i * steps) : ((total * i + tden) * steps)] , args=(Tden, dGs))

        concentration[(total * i * steps) : ((total * i + tden) * steps)] = integration_den


        # dGs[0] = amplicon_length * dH_S1S2 - Tanneal * dS
        #
        # dGs[1] = primer_length * dH_S1P2_S2P1 - Tanneal * dS
        #
        # dGs[2] = extended_primer * dH_S1Q2_S2Q1 - Tanneal * dS
        #
        # dGs[3] = dH_enzyme - Tanneal * dS



        dGs[0] = amplicon_length * dH - Tanneal * dS

        dGs[1] = primer_length * dH - Tanneal * dS

        dGs[2] = extended_primer * dH - Tanneal * dS

        dGs[3] = dH_enzyme2 - Tanneal * dS

        dGs = np.clip(dGs, a_min= None, a_max= 1e+14 )

        print("dGs_den", dGs)




        print("dGs_anneals", dGs)

        integration_anneal = odeint(PCR_reaction, integration_den[-1], time[((total * i + tden) * steps)-1 : ((total * i + tden + tanneal) * steps)] , args=(Tanneal, dGs ))

        concentration[((total * i + tden) * steps) - 1 : ((total * i + tden + tanneal) * steps)] = integration_anneal


        # dGs[0] = amplicon_length * dH_S1S2 - Text * dS
        #
        # dGs[1] = primer_length * dH_S1P2_S2P1 - Text * dS
        #
        # dGs[2] = extended_primer * dH_S1Q2_S2Q1 - Text * dS
        #
        # dGs[3] = dH_enzyme - Text * dS



        dGs[0] = amplicon_length * dH - Text * dS

        dGs[1] = primer_length * dH - Text * dS

        dGs[2] = extended_primer * dH - Text * dS

        dGs[3] = dH_enzyme2 - Text * dS

        dGs = np.clip(dGs, a_min= None, a_max= 1e+14 )




        print("dGs_text", dGs)

        integration_ext = odeint(PCR_reaction, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) -1 : (total * (i + 1) * steps)], args=(Text, dGs))

        concentration[((total * i + tden + tanneal) * steps) -1 : (total * (i + 1) * steps)] = integration_ext


        values = integration_ext[-1]



    print("The concentration of the 17 species at the end of the process:", values)

    print("concentration", concentration)


    #Plotting the concentrartions over time in two figures

    plt.figure(1)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)

    plots = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]


    for i in range(10):

        plt.subplot(2, 5, i+1)


        plt.plot(time, concentration[:, plots[i]])


        plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time")
        plt.ylabel("Concentration")


    plt.figure(2)

    for i in range(8):

        plt.subplot(2, 4, i+1)


        plt.plot(time, concentration[:, i])

        plt.legend([species[i]], loc='best', prop={'size':10})

        plt.xlabel("Time")
        plt.ylabel("Concentration")


    plt.figure(3)

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


    # overflow encountered in exp
    # kr2a = kf2 * np.exp(dGs[1]/(R*T))
    # RuntimeWarning: invalid value encountered in double_scalars
    # rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2
