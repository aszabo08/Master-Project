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



bases = ["A", "G", "C", "U"]

amplicon_list = []

for i in range(1000):

    amplicon_list.append(random.choice(bases))


amplicon = "".join(amplicon_list)


print(amplicon)


primer_list = []

for i in range(15):

    primer_list.append(random.choice(bases))


primer = "".join(primer_list)


print(primer)






"The species included in the reactions are the following: "

values = [S1S2, S1, S2, P1, P2, S1P2, S2P1, E, S1P2E, S2P1E, dNTP, S1P2EdNTP, S2P1EdNTP]


kB = 1.38064852 * (1/np.power(10,23))


def denaturation(Tden):

    "used sequences: S1:AAGGCCCUAAU  complementary S2: UUCCGGGAUUA for RNAcofold"


#    kf1 = 1

#    dG = -1316.7048

#   kB = 1.38064852 * (1/np.power(10,23))

#   kr1 = kf1 * np.exp(dG/kB*Tden)



# or: S1S2/ S1 * S2 = kf1 / kr1

#   kr1 = (kf1 * S1 * S2 )/ S1S2

    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]


    kf1 = 1

    dG = (np.log(S1S2) - np.log(S1 * S2)) * kB * Tden

    kr1 = kf1 * np.exp(dG/kB*Tden)


    dS1S2dt = -kf1 * S1S2 + kr1 * S1 * S2
    dS1dt = kf1 * S1S2 - kr1 * S1 * S2
    dS2dt = kf1 * S1S2 - kr1 * S1 * S2

    y = np.empty(3)

    y[0] = dS1S2dt
    y[1] = dS1dt
    y[2] = dS2dt

    return y


def primer_binding(Tanneal):


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

    # dG is equal for the two reactions as S1 = S2 and P1 = p2?

    dGa = (np.log(S1P2) - np.log(S1 * P2)) * kB * Tanneal

    kr2a = kf2 * np.exp(dGa/kB*Tanneal)


    dS1dt = -kf2 * S1 * P2 + kr2a * S1P2
    dP2dt = -kf2 * S1 * P2 + kr2a * S1P2
    dS1P2dt = kf2 * S1 * P2 - kr2a * S1P2


    dGb = (np.log(S2P1) - np.log(S2 * P1)) * kB * Tanneal

    kr2b = kf2 * np.exp(dGb/kB*Tanneal)


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



def ploymerase_binding(T):


    SP = values[0]
    E = values[1]
    SPE = values[2]


    dSPdt = -kf3 * SP * E + kr3 * SPE

    dEdt = -kf3 * SP * E + kr3 * SPE

    dSPEdt = kf3 * SP * E - kr3 * SPE

    y = np.empty(3)

    y[0] = dSPdt
    y[1] = dEdt
    y[2] = dSPEdt

    return y


def stabilizing(values, t, kf4, kr4):


    SPE = values[0]
    N = values[1]
    SPEN = values[2]

    dSPEdt = - kf4 * SPE * N + kr4 * SPEN

    dNdt = - kf4 * SPE * N + kr4 * SPEN

    dSPENdt = kf4 * SPE * N - kr4 * SPEN

    y[0] = dSPEdt
    y[1] = dNdt
    y[2] = dSPENdt

    return y


def primer_extension(values, t, ke):

    "N_comp = complementary nucleotides in the whole length of the substrate"

    SPEN = values[0]
    N_comp = values[1]


    dSPNN_compdt + dEdt = ke * SPEN * N_comp

    y = dSPNN_compdt + dEdt

    return y
















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



