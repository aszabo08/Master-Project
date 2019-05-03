"""
This application visualizes primer binding and misbinding to a substrate, creating substrate-primer products.
"S" stands for the concentration of the substrate, "P" and "SP" for the concentration of the primer and the bound state of the substrate and primer.
The concentrations for all the above mentioned species are given in mole/liter.
The k parameters (k1, k2, k3) are rate constants in binding and misbinding processes.
"""


from scipy.integrate import odeint

import matplotlib.pyplot as plt

import numpy as np

import random

amplicon_length = 1000

primer_length = 15

extended_length = amplicon_length - primer_length







def denaturation(values, t, Tden):


    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]


    kf1 = 1
    kr1 = kf1 * np.exp(dG_den/kB*Tden)


    dS1S2dt = -kf1 * S1S2 + kr1 * S1 * S2
    dS1dt = kf1 * S1S2 - kr1 * S1 * S2
    dS2dt = kf1 * S1S2 - kr1 * S1 * S2

    y = np.empty(3)

    y[0] = dS1S2dt
    y[1] = dS1dt
    y[2] = dS2dt

    return y



def primer_binding(values, t, Tanneal):


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
    kr2a = kf2 * np.exp(dG_bind1/kB*Tanneal)


    dS1dt = -kf2 * S1 * P2 + kr2a * S1P2
    dP2dt = -kf2 * S1 * P2 + kr2a * S1P2
    dS1P2dt = kf2 * S1 * P2 - kr2a * S1P2


    kr2b = kf2 * np.exp(dG_bind2/kB*Tanneal)


    dS2dt = -kf2 * S2 * P1 + kr2b * S2P1
    dP1dt = -kf2 * S2 * P1 + kr2b * S2P1
    dS2P1dt = kf2 * S2 * P1 - kr2b * S2P1

    y = np.empty(6)

    y[0] = dS1dt
    y[1] = dS2dt
    y[2] = dP1dt
    y[3] = dP2dt
    y[4] = dS1P2dt
    y[5] = dS2P1dt

    return y


def ploymerase_binding(values, t, Text):


    S1P2 = values[5]
    S2P1 = values[6]
    E = values[7]
    S1P2E = values[8]
    S2P1E = values[9]


    kf3 = 1


    kr3a = kf3 * np.exp(dGa/kB*Text)

    dS1P2dt = -kf3 * S1P2 * E + kr3a * S1P2E

    dEadt = -kf3 * S1P2 * E + kr3a * S1P2E

    dS1P2Edt = kf3 * S1P2 * E - kr3a * S1P2E



    kr3b = kf3 * np.exp(dGb/kB*Text)

    dS2P1dt = -kf3 * S2P1 * E + kr3b * S2P1E

    dEbdt = -kf3 * S2P1 * E + kr3b * S2P1E

    dS2P1Edt = kf3 * S2P1 * E - kr3b * S2P1E



    y = np.empty(6)

    y[0] = dS1P2dt
    y[1] = dS2P1dt
    y[2] = dEadt
    y[3] = dEbdt
    y[4] = dS1P2Edt
    y[5] = dS2P1Edt


    return y


def stabilizing(values, t, T):

    """ The primer is extended by a few nucleotides to ensure the primer binding
        to the substrate without dissociation when reaching the extension temperature

        In this model 99.5 % of the primer - substrate complexes stay together, while 0.5 % of them will melt
        at the extension temperature"""


    S1P2E = values[8]
    S2P1E = values[9]
    dNTP = values[10]
    Q1 = values[11]         # Q1 = P1 + n * dNTP
    Q2 = values[12]         # Q2 = P2 + n * dNTP
    S1Q2E = values[13]
    S2Q1E = values[14]
    


   # 0.995 = np.exp(-(dG/kB*T))

    dG = - np.log(0.995) * kB * Text

    dG = n * dH - T * dS

    n = 5       # number of bases in the initial extension


    #   S1P2E + n * dNTP -> S1Q2E

    ce = 1000   # [dNTP/s] concentration of polymerase enzyme

    #   s1 = ce / n * S1P2E * dNTP      # the speed of the reaction

    dS1P2Edt = -(ce / n * S1P2E * dNTP)

    dS1Q2Edt = ce / n * S1P2E * dNTP

    #   S2P1E + n * dNTP -> S2Q1E

    #   s2 = ce / n * S2P1E * dNTP      # the speed of the reaction

    dS2P1Edt = -(ce / n * S2P1E * dNTP )

    dS2Q1Edt = ce / n * S2P1E * dNTP


    ddNTPdt =  -(ce * S1P2E * dNTP) -(ce * S2P1E * dNTP )       # -n * (ce / n * S1P2E * dNTP) and -n * (ce / n * S2P1E * dNTP )


    y = np.empty(5)

    y[0] = dS1P2Edt
    y[1] = dS2P1Edt
    y[2] = ddNTPdt
    y[3] = dS1Q2Edt
    y[4] = dS2Q1Edt


    return y


def primer_extension(values, t, Text):


    extended_primer_length = n + len(primer)

    extension = len(amplicon) - extended_primer_length


    a = (ke * S1P2E * dNTP) / extension


    return





if __name__ == '__main__':



    dG_den = -2573.5784                             # using RNAcofold with the two complementary sequences: dG = -615.10 kcal/mol which is equivalent with -2573.5784 kJ/mol

    dG_bind1 = -1306.6632                           # using RNAcofold with amplicon and primer: dG = -312.30 kcal/mol which is equivalent with -1306.6632 kJ/mol

    dG_bind2 = -1294.948                                  # using RNAcofold with the complement of the amplicon and the complement of the primer: dG = -309.50 kcal/mol which is equivalent with -1294.948 kJ/mol

    kB = 1.38064852 * (1/np.power(10, 23))          # Boltzmann constant"


    # "The species included in the reactions are the following: "

    # state vector:

    # values = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2, "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]








    for i in range(number_cycles):


        t1 = np.linspace(0, tden, tden * 10)

        int_denaturation = odeint(denaturation, values[0:3], t1, args=(Tden, ))

        values[1] = list(int_denaturation[-1])[1]

        values[2] = list(int_denaturation[-1])[2]

        t2 = np.linspace(tden, tden + tanneal, tanneal * 10)

        int_binding = odeint(primer_binding, values[1:7], t2, args=(Tanneal, ))











def primer_misbinding(values, t, kfm, kram, krbm):

    """
    Describes the differential equations for S, P, SP1 and SP2 in a dynamic biological reaction

    Args:
        values: value of S, P, SP1 (bound state of S and P, with correct binding location)
                and SP2 (bound state of S and P, with incorrect binding location) after t time period
        t: time period
        k1: rate constant for binding ( S + P -> SP1, S + P -> SP2)
        k2: rate constant for dissociation, when the binding location is correct ( SP1 -> S + P)
        k3: rate constant for dissociation, when the binding location is incorrect ( SP2 -> S + P)

    Returns:
        An array with the result of the differential equations for S, P, SP1 and SP2 over t time
    """

    S = values[0]
    P = values[1]
    SP1 = values[2]
    SP2 = values[3]


    dSdt = -2 * kfm * S * P + kram * SP1 + krbm * SP2
    dPdt = -2 * kfm * S * P + kram * SP1 + krbm * SP2
    dSP1dt = kfm * S * P - kram * SP1
    dSP2dt = kfm* S * P - krbm * SP2

    y = np.empty(4)

    y[0] = dSdt
    y[1] = dPdt
    y[2] = dSP1dt
    y[3] = dSP2dt

    return y



