

import PCR_project

new species: S1M2

extend values

change y size


# amplicon length - primer length - n = 1000 - 15 -10 = 975         -------->  the longest extended length is 975 -1, when the product is only one nucleotide short
#                                                                              the shortest extended length is 0
# average extended length = (974 + 0)/2 = 487

misbinding_extended_length = 487





def misbinding_primer_1(values, t, T, dGs):



    S1 = values[1]
    P2 = values[4]
    S1M2 = values[17]


    kf6 = forward_rate

    #exponent_2a = dGs[1]/(R*T)

    #exponent_2a = np.clip((dGs[1]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_6a = exponent_clipping(dGs[1]/(R*T))

    kr6a = kf6 * np.exp(exponent_6a)

    kf6, kr6a = clipping(kf6, kr6a)

    # if np.abs(kr2a < 1e-14):
    #     kr2a = 0

    #rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2

    #rate_S1P2_bind = np.clip((kf2 * S1 * P2 - kr2a * S1P2), a_min= min_clip, a_max= max_clip)

    rate_S1M2_bind = rate_clipping(kf6 * S1 * P2 - kr6a * S1M2)


    y = np.zeros(17)


    y[1] = - rate_S1M2_bind
    y[4] = - rate_S1M2_bind
    y[17] = rate_S1M2_bind



    return y





def misbinding_polymerase_1(values, t, T, dGs):


    S1M2 = values[17]
    E = values[7]
    S1M2E = values[18]


    kf7 = forward_rate

    #exponent_3a = dGs[3]/(R*T)

    #exponent_3a = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_7a = exponent_clipping(dGs[3]/(R*T))

    kr7a = kf7 * np.exp(exponent_7a)

    kf7, kr7a = clipping(kf7, kr7a)

    #rate_poly_S1P2_bind = kf3 * S1P2 * E - kr3a * S1P2E

    #rate_poly_S1P2_bind = np.clip((kf3 * S1P2 * E - kr3a * S1P2E), a_min= min_clip, a_max= max_clip)

    rate_poly_S1M2_bind = rate_clipping(kf7 * S1M2 * E - kr7a * S1M2E)

    y = np.zeros(17)

    y[17] = - rate_poly_S1M2_bind
    y[7] = - rate_poly_S1M2_bind
    y[18] = rate_poly_S1M2_bind



    return y





def misbinding_primer_ext_1(values, t, T, dGs):

    """ The primer is extended by a few nucleotides to ensure the primer binding
        to the substrate without dissociation when reaching the extension temperature

        In this model 99.5 % of the primer - substrate complexes stay together, while 0.5 % of them will melt
        at the extension temperature"""


    S1M2E = values[18]
    dNTP = values[10]

    #if (total_molecule_length(dNTP, 1)) >= (2 * n):


        #primer_ext_1.counter += 1


    #   S1P2E + n * dNTP -> S1Q2E
    #   S2P1E + n * dNTP -> S2Q1E

    #ce = 1000   # [dNTP/s] concentration of polymerase enzyme

    ce = taq_nt_per_s(T)

    #rate_ext_1 = (ce / n) * S1P2E * dNTP

    #rate_ext_1 = np.clip(((ce / n) * S1P2E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_1 = rate_clipping((ce / n) * S1M2E * dNTP)

    y = np.zeros(17)

    y[18] = - rate_ext_1  # concentration of S1M2E
    y[10] = - n * rate_ext_1  # concentration of dNTP
    y[19] = rate_ext_1  # concentration of S1N2E


    return y




def misbinding_primer_ext_2(values, t, T, dGs):



    S1N2E = values[19]
    dNTP = values[10]

    ce_N = taq_nt_per_s(T)

    # reaction: S1N2E + misbinding_extended_length * dNTP ---> L1L2 + E

    #rate_ext_Q1 = (ce_Q / extended_length) * S1Q2E * dNTP

    #rate_ext_Q1 = np.clip(((ce_Q / extended_length) * S1Q2E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_N1 = rate_clipping((ce_N / misbinding_extended_length) * S1N2E * dNTP)


    y = np.zeros(17)


    y[7] = rate_ext_N1
    y[10] = - misbinding_extended_length * rate_ext_N1           # concentration of dNTP
    y[19] = - rate_ext_N1       # concentration of S1N2E
    y[20] = rate_ext_N1

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






def misbinding_denaturation(values, t, T, dGs):


    S1L2 = values[20]
    S1 = values[1]
    L2 = values[21]


    kf8 = forward_rate


    #exponent_1 = dGs[0]/(R*T)

    exponent_8 = exponent_clipping(dGs[0]/(R*T))

    kr8 = kf8 * np.exp(exponent_8)


    kf8, kr8 = clipping(kf8, kr8)



    #rate_den = - kr1 * S1S2 + kf1 * S1 * S2

    #rate_den = np.clip((- kr1 * S1S2 + kf1 * S1 * S2), a_min= min_clip, a_max= max_clip)

    rate_den = rate_clipping(- kr8 * S1L2 + kf8 * S1 * L2)


    y = np.zeros(17)

    y[20] = rate_den

    y[1] = -rate_den

    y[21] = -rate_den



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
