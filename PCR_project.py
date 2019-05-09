"""
This application visualizes primer binding and misbinding to a substrate, creating substrate-primer products.
"S" stands for the concentration of the substrate, "P" and "SP" for the concentration of the primer and the bound state of the substrate and primer.
The concentrations for all the above mentioned species are given in mole/liter.
The k parameters (k1, k2, k3) are rate constants in binding and misbinding processes.
"""


from scipy.integrate import odeint

#import matplotlib.pyplot as plt

import numpy as np

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

    y = np.empty(17)

    y[0] = dS1S2dt
    y[1] = dS1dt
    y[2] = dS2dt
    y[3] = values[3]    # concentration of P1 in microL
    y[4] = values[4]    # concentration of P2 in microL
    y[5] = values[5]      # concentration of S1P2
    y[6] = values[6]     # concentration of S2P1
    y[7] = values[7]     # concentration of E
    y[8] = values[8]      # concentration of S1P2E
    y[9] = values[9]     # concentration of S2P1E
    y[10] = values[10]    # concentration of dNTP
    y[11] = values[11]      # concentration of Q1
    y[12] = values[12]     # concentration of Q2
    y[13] = values[13]       # concentration of S1Q2E
    y[14] = values[14]      # concentration of S2Q1E
    y[15] = values[15]       #concentration of S2Q1
    y[16] = values[16]       #concentration of S1Q2

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


    y = np.empty(17)


    y[0] = values[0]
    y[1] = dS1dt
    y[2] = dS2dt
    y[3] = dP1dt
    y[4] = dP2dt
    y[5] = dS1P2dt
    y[6] = dS2P1dt
    y[7] = values[7]     # concentration of E
    y[8] = values[8]      # concentration of S1P2E
    y[9] = values[9]     # concentration of S2P1E
    y[10] = values[10]    # concentration of dNTP
    y[11] = values[11]      # concentration of Q1
    y[12] = values[12]     # concentration of Q2
    y[13] = values[13]       # concentration of S1Q2E
    y[14] = values[14]      # concentration of S2Q1E
    y[15] = values[15]       #concentration of S2Q1
    y[16] = values[16]      #concentration of S1Q2


    return y


def polymerase_binding(values, t, Text):


    S1P2 = values[5]
    S2P1 = values[6]
    E = values[7]
    S1P2E = values[8]
    S2P1E = values[9]


    kf3 = 1


    kr3a = kf3 * np.exp(dGa/kB*Text)

    dS1P2dt = -kf3 * S1P2 * E + kr3a * S1P2E

    dS1P2Edt = kf3 * S1P2 * E - kr3a * S1P2E


    kr3b = kf3 * np.exp(dGb/kB*Text)

    dS2P1dt = -kf3 * S2P1 * E + kr3b * S2P1E

    dS2P1Edt = kf3 * S2P1 * E - kr3b * S2P1E

    dEdt = (-kf3 * S1P2 * E + kr3a * S1P2E) + (-kf3 * S2P1 * E + kr3b * S2P1E)

    # dEadt = (-kf3 * S1P2 * E + kr3a * S1P2E)

    # dEbdt = -kf3 * S2P1 * E + kr3b * S2P1E


    y = np.empty(17)


    y[0] = values[0]
    y[1] = values[1]
    y[2] = values[2]
    y[3] = values[3]
    y[4] = values[4]
    y[5] = dS1P2dt
    y[6] = dS2P1dt
    y[7] = dEdt
    y[8] = dS1P2Edt
    y[9] = dS2P1Edt
    y[10] = values[10]    # concentration of dNTP
    y[11] = values[11]      # concentration of Q1
    y[12] = values[12]     # concentration of Q2
    y[13] = values[13]       # concentration of S1Q2E
    y[14] = values[14]      # concentration of S2Q1E
    y[15] = values[15]       #concentration of S2Q1
    y[16] = values[16]      #concentration of S1Q2



    return y


def stabilizing(values, t, T):

    """ The primer is extended by a few nucleotides to ensure the primer binding
        to the substrate without dissociation when reaching the extension temperature

        In this model 99.5 % of the primer - substrate complexes stay together, while 0.5 % of them will melt
        at the extension temperature"""


    S1P2E = values[8]
    S2P1E = values[9]
    dNTP = values[10]
    


    # 0.995 = np.exp(-(dG/kB*T))

    dG = - np.log(0.995) * kB * Text

    dG = n * dH - T * dS

    n = 5       # number of bases in the initial extension


    #   S1P2E + n * dNTP -> S1Q2E

    ce = 1000   # [dNTP/s] concentration of polymerase enzyme

    #   s_1 = ce / n * S1P2E * dNTP      # the speed of the reaction creating complexes

    dS1P2Edt = -(ce / n * S1P2E * dNTP) - 3 * (ce2 / n * S1P2E * dNTP)

    dS1Q2Edt = ce / n * S1P2E * dNTP


    # ce2 is the speed of creating S1, G2 and E which are going to be denaturated

    # s_1 = ce2 / n * S1P2E * dNTP

    dS1dt = ce2 / n * S1P2E * dNTP

    dQ2dt = ce2 / n * S1P2E * dNTP

    # dEdt = ce2 / n * S1P2E * dNTP




    #   S2P1E + n * dNTP -> S2Q1E

    #   s_ 2 = ce2b / n * S2P1E * dNTP      # the speed of the reaction creating complexes

    dS2P1Edt = -(ce / n * S2P1E * dNTP ) -3 * (ce2 / n * S1P2E * dNTP)

    dS2Q1Edt = ce / n * S2P1E * dNTP

    # s_2a = ce2 / n * S2P1E * dNTP

    dS2dt = ce2 / n * S2P1E * dNTP

    dQ1dt = ce2 / n * S2P1E * dNTP

    dEdt = ce2 / n * S1P2E * dNTP + ce2 / n * S2P1E * dNTP


    ddNTPdt =  -(ce * S1P2E * dNTP) -(ce * S2P1E * dNTP )- 3 * (ce2 / n * S1P2E * dNTP) -3 * (ce2 / n * S1P2E * dNTP)

    y = np.empty(17)

    y[0] = values[0]
    y[1] = dS1dt
    y[2] = dS2dt
    y[3] = dP1dt
    y[4] = dP2dt
    y[4] = dS1P2dt
    y[5] = dS2P1dt
    y[6] = values[6]    # concentration of S2P1
    y[7] = dEdt         # concentration of E
    y[8] = dS1P2Edt     # concentration of S1P2E
    y[9] = dS2P1Edt     # concentration of S2P1E
    y[10] = ddNTPdt     # concentration of dNTP
    y[11] = dQ1dt       # concentration of Q1
    y[12] = dQ2dt       # concentration of Q2
    y[13] = dS1Q2Edt    # concentration of S1Q2E
    y[14] = dS2Q1Edt    # concentration of S2Q1E
    y[15] = values[15]       #concentration of S2Q1
    y[16] = values[16]      #concentration of S1Q2


    return y


def primer_extension(values, t, Text):


    extended_primer_length = n + len(primer)

    extension = len(amplicon) - extended_primer_length


    a = (ke * S1P2E * dNTP) / extension


    return





if __name__ == '__main__':


    values = [0 for i in range(17)]

    values[0] = 60      # concentration of plasmid (S1S2) in ng
    values[1] = 0       # concentration of S1 in ng
    values[2] = 0       # concentration of S2 in ng
    values[3] = 9.2     # concentration of P1 in microL
    values[4] = 9.2     # concentration of P2 in microL
    values[5] = 0       # concentration of S1P2
    values[6] = 0       # concentration of S2P1
    values[7] = 8.2      # concentration of E
    values[8] = 0       # concentration of S1P2E
    values[9] = 0       # concentration of S2P1E
    values[10] = 4.2    # concentration of dNTP
    values[11] = 0       # concentration of Q1
    values[12] = 0       # concentration of Q2
    values[13] = 0       # concentration of S1Q2E
    values[14] = 0       # concentration of S2Q1E
    values[15] = 0       #concentration of S2Q1
    values[16] = 0       #concentration of S1Q2



    Tden = 90           # degree Celsius
    Tanneal = 72        # degree Celsius
    Text = 80           # degree Celsius

    tden = 2           # seconds
    tanneal = 3        # seconds
    text = 4           # seconds

    number_cycles = 20


    time = np.linspace(0, tden + tanneal + text, (tden + tanneal + text) * 10)      # for every second 10 points are distinguished

    #  time = np.array([(np.linspace(0, tden, tden * 10)), (np.linspace(tden, tden + tanneal, tanneal * 10)), (np.linspace(tden + tanneal, tden + tanneal + text, text * 10))])

    T = np.array([([Tden]*tden*10), ([Tanneal]*tanneal*10), ([Text]*text*10)])


    dG_den = -2573.5784                             # using RNAcofold with the two complementary sequences: dG = -615.10 kcal/mol which is equivalent with -2573.5784 kJ/mol

    dG_bind1 = -1306.6632                           # using RNAcofold with amplicon and primer: dG = -312.30 kcal/mol which is equivalent with -1306.6632 kJ/mol

    dG_bind2 = -1294.948                                  # using RNAcofold with the complement of the amplicon and the complement of the primer: dG = -309.50 kcal/mol which is equivalent with -1294.948 kJ/mol

    kB = 1.38064852 * (1/np.power(10, 23))          # Boltzmann constant"


    # "The species included in the reactions are the following: "

    # state vector:

    # values = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2, "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]

    change = np.zeros((number_cycles, 17), dtype=np.int)


    for i in range(number_cycles):


        change[i] = odeint(denaturation, values, time, args=(T, ))[-1] + odeint(primer_binding, values, time, args=(T, ))[-1] + odeint(polymerase_binding, values, time, args=(T, ))[-1] + odeint(stabilizing, values, time, args=(T, ))[-1]

        values = values + change[i]

result = values



















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


        #
        # values[1] = list(int_denaturation[-1])[1]
        #
        # values[2] = list(int_denaturation[-1])[2]
