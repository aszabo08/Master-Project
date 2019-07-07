

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
    y[20] = rate_ext_N1         # concentration of S1L2

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

    y[20] = rate_den        # concentration of S1L2

    y[1] = -rate_den        # concentration of S1

    y[21] = -rate_den       # concentration of L2



    return y



def misbinding_polymerase_2(values, t, T, dGs):



    E = values[7]
    S1N2E = values[19]
    S1N2 = values[22]



    kf9 = forward_rate

    #exponent_4a = dGs[3]/(R*T)

    #exponent_4a = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_9a = exponent_clipping(dGs[3]/(R*T))


    kr9a = kf9 * np.exp(exponent_9a)

    kf9, kr9a = clipping(kf9, kr9a)

    #rate_poly_S1Q2_bind = kf4 * S1Q2 * E - kr4a * S1Q2E

    #rate_poly_S1Q2_bind = np.clip((kf4 * S1Q2 * E - kr4a * S1Q2E), a_min= min_clip, a_max= max_clip)

    rate_poly_S1N2_bind = rate_clipping(kf9 * S1N2 * E - kr9a * S1N2E)

    y = np.zeros(17)

    y[22] = - rate_poly_S1N2_bind
    y[7] = - rate_poly_S1N2_bind
    y[19] = rate_poly_S1N2_bind



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






def misbinding_primer_2(values, t, T, dGs):


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
    N2 = values[23]
    S1N2 = values[22]



    kf10 = forward_rate

    #exponent_5a = dGs[2]/(R*T)

    #exponent_5a = np.clip((dGs[2]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_10a = exponent_clipping(dGs[2]/(R*T))

    kr10a = kf10 * np.exp(exponent_10a)

    kf10, kr10a = clipping(kf10, kr10a)

    #rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    rate_S1N2_bind = rate_clipping(kf10 * S1 * N2 - kr10a * S1N2)

    y = np.zeros(17)


    y[1] = - rate_S1N2_bind
    y[23] = - rate_S1N2_bind
    y[22] = rate_S1N2_bind



    return y















def short_misbinding_primer(values, t, T, dGs):


    L2P1 = values[24]
    P1 = values[3]
    L2 = values[21]


    L1P2 = values[30]
    P2 = values[4]
    L1 = values[29]



    kf11 = forward_rate


    #exponent_1 = dGs[0]/(R*T)

    exponent_11a = exponent_clipping(dGs[0]/(R*T))

    kr11a = kf11 * np.exp(exponent_11a)


    kf11, kr11a = clipping(kf11, kr11a)

    rate_L2P1 = rate_clipping(- kr11a * L2P1 + kf11 * P1 * L2)



    #exponent_1 = dGs[0]/(R*T)

    exponent_11b = exponent_clipping(dGs[0]/(R*T))

    kr11b = kf11 * np.exp(exponent_11b)


    kf11, kr11b = clipping(kf11, kr11b)

    rate_L1P2 = rate_clipping(- kr11b * L1P2 + kf11 * P2 * L1)



    y = np.zeros(17)

    y[24] = rate_L2P1        # concentration of P1L2

    y[3] = -rate_L2P1        # concentration of P1

    y[21] = -rate_L2P1       # concentration of L2





    y[30] = rate_L1P2        # concentration of L1P2

    y[4] = -rate_L1P2        # concentration of P2

    y[29] = -rate_L1P2      # concentration of L1



    return y




def short_misbinding_polymerase_1(values, t, T, dGs):


    L2P1 = values[24]
    E = values[7]
    L2P1E = values[25]



    L1P2 = values[30]
    #E = values[7]
    L1P2E = values[31]



    kf12 = forward_rate

    #exponent_3a = dGs[3]/(R*T)

    #exponent_3a = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_12a = exponent_clipping(dGs[3]/(R*T))

    kr12a = kf12 * np.exp(exponent_12a)

    kf12, kr12a = clipping(kf12, kr12a)

    rate_poly_L2P1_bind = rate_clipping(kf12 * L2P1 * E - kr12a * L2P1E)






    #exponent_3a = dGs[3]/(R*T)

    #exponent_3a = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_12b = exponent_clipping(dGs[3]/(R*T))

    kr12b = kf12 * np.exp(exponent_12b)

    kf12, kr12b = clipping(kf12, kr12b)

    rate_poly_L1P2_bind = rate_clipping(kf12 * L1P2 * E - kr12b * L1P2E)


    enzyme_misbinding = rate_clipping(- rate_poly_L2P1_bind - rate_poly_L1P2_bind)


    y = np.zeros(17)

    y[24] = - rate_poly_L2P1_bind
    #y[7] = - rate_poly_L2P1_bind
    y[25] = rate_poly_L2P1_bind


    y[30] = - rate_poly_L1P2_bind
    y[7] = enzyme_misbinding
    y[31] = rate_poly_L1P2_bind



    return y



def short_misbinding_primer_ext_1(values, t, T, dGs):

    """ The primer is extended by a few nucleotides to ensure the primer binding
        to the substrate without dissociation when reaching the extension temperature

        In this model 99.5 % of the primer - substrate complexes stay together, while 0.5 % of them will melt
        at the extension temperature"""


    L2P1E = values[25]
    dNTP = values[10]


    L1P2E = values[31]




    #if (total_molecule_length(dNTP, 1)) >= (2 * n):


        #primer_ext_1.counter += 1


    #   S1P2E + n * dNTP -> S1Q2E
    #   S2P1E + n * dNTP -> S2Q1E

    #ce = 1000   # [dNTP/s] concentration of polymerase enzyme

    ce = taq_nt_per_s(T)

    #rate_ext_1 = (ce / n) * S1P2E * dNTP

    #rate_ext_1 = np.clip(((ce / n) * S1P2E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_1 = rate_clipping((ce / n) * L2P1E * dNTP)


    rate_ext_2 = rate_clipping((ce / n) * L1P2E * dNTP)


    misbinding_nucleotide = rate_clipping(- n * rate_ext_1 - n * rate_ext_2)

    y = np.zeros(17)

    y[25] = - rate_ext_1                 # concentration of L2P1E
    #y[10] = - n * rate_ext_1             # concentration of dNTP
    y[26] = rate_ext_1                   # concentration of L2Q1E


    y[31] = - rate_ext_2                 # concentration of L1P2E
    y[10] = misbinding_nucleotide            # concentration of dNTP
    y[32] = rate_ext_2                   # concentration of L1Q2E


    return y




def short_misbinding_primer_ext_2(values, t, T, dGs):



    L2Q1E = values[26]
    dNTP = values[10]


    L1Q2E = values[32]



    ce_L = taq_nt_per_s(T)

    # reaction: S1N2E + misbinding_extended_length * dNTP ---> L1L2 + E

    #rate_ext_Q1 = (ce_Q / extended_length) * S1Q2E * dNTP

    #rate_ext_Q1 = np.clip(((ce_Q / extended_length) * S1Q2E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_QL_1 = rate_clipping((ce_L / misbinding_extended_length) * L2Q1E * dNTP)




    rate_ext_QL_2 = rate_clipping((ce_L / misbinding_extended_length) * L1Q2E * dNTP)


    mis_product = rate_clipping(rate_ext_QL_1 + rate_ext_QL_2)

    mis_nucleotide = rate_clipping(- misbinding_extended_length * rate_ext_QL_1- misbinding_extended_length * rate_ext_QL_2)


    y = np.zeros(17)


    #y[7] = rate_ext_QL_1          # concentration of enzyme
    #y[10] = - misbinding_extended_length * rate_ext_QL_1           # concentration of dNTP
    y[26] = - rate_ext_QL_1       # concentration of L2Q1E
    #y[27] = rate_ext_QL_1         # concentration of L1L2



    y[7] = mis_product          # concentration of enzyme
    y[10] = mis_nucleotide          # concentration of dNTP
    y[32] = - rate_ext_QL_2      # concentration of L1Q2E
    y[27] = mis_product        # concentration of L1L2


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





def short_misbinding_polymerase_2(values, t, T, dGs):



    E = values[7]
    L2Q1E = values[26]
    L2Q1 = values[28]


    L1Q2E = values[32]
    L1Q2 = values[33]



    kf13 = forward_rate

    #exponent_4a = dGs[3]/(R*T)

    #exponent_4a = np.clip((dGs[3]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_13a = exponent_clipping(dGs[3]/(R*T))


    kr13a = kf13 * np.exp(exponent_13a)

    kf13, kr13a = clipping(kf13, kr13a)

    rate_poly_Q1L2_bind = rate_clipping(kf13 * L2Q1 * E - kr13a * L2Q1E)





    exponent_13b = exponent_clipping(dGs[3]/(R*T))


    kr13b = kf13 * np.exp(exponent_13b)

    kf13, kr13b = clipping(kf13, kr13b)

    rate_poly_L1Q2_bind = rate_clipping(kf13 * L1Q2 * E - kr13b * L1Q2E)



    y = np.zeros(17)

    y[28] = - rate_poly_Q1L2_bind
    y[7] = - rate_poly_Q1L2_bind
    y[26] = rate_poly_Q1L2_bind



    y[33] = - rate_poly_L1Q2_bind
    y[7] = - rate_poly_L1Q2_bind
    y[32] = rate_poly_L1Q2_bind




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






def short_misbinding_primer_2(values, t, T, dGs):


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


    Q1 = values[11]
    L2 = values[21]
    L2Q1 = values[28]



    Q2 = values[12]
    L1 = values[29]
    L1Q2 = values[33]



    kf14 = forward_rate

    #exponent_5a = dGs[2]/(R*T)

    #exponent_5a = np.clip((dGs[2]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_14a = exponent_clipping(dGs[2]/(R*T))

    kr14a = kf14 * np.exp(exponent_14a)

    kf14, kr14a = clipping(kf14, kr14a)

    #rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    rate_Q1L2_bind = rate_clipping(kf14 * Q1 * L2 - kr14a * L2Q1)





    #exponent_5a = dGs[2]/(R*T)

    #exponent_5a = np.clip((dGs[2]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_14b = exponent_clipping(dGs[2]/(R*T))

    kr14b = kf14 * np.exp(exponent_14b)

    kf14, kr14b = clipping(kf14, kr14b)

    #rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    rate_L1Q2_bind = rate_clipping(kf14 * Q2 * L1 - kr14b * L1Q2)




    y = np.zeros(17)


    y[11] = - rate_Q1L2_bind
    y[21] = - rate_Q1L2_bind
    y[28] = rate_Q1L2_bind





    y[12] = - rate_L1Q2_bind
    y[29] = - rate_L1Q2_bind
    y[33] = rate_L1Q2_bind





    return y






def L_misbinding_denaturation(values, t, T, dGs):


    L1L2 = values[27]
    L1 = values[29]
    L2 = values[21]


    kf15 = forward_rate


    #exponent_1 = dGs[0]/(R*T)

    exponent_15 = exponent_clipping(dGs[0]/(R*T))

    kr15 = kf15 * np.exp(exponent_15)


    kf15, kr15 = clipping(kf15, kr15)

    rate_den = rate_clipping(- kr15 * L1L2 + kf15 * L1 * L2)


    y = np.zeros(17)

    y[27] = rate_den        # concentration of L1L2

    y[29] = -rate_den        # concentration of L1

    y[21] = -rate_den       # concentration of L2



    return y







def PCR_reaction(values, t, T, dGs):  # not using rate clipping cos the array data

    summary = denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + primer_binding_2(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + primer_ext_2(values, t, T, dGs) + taq_denaturation(values, t, T, dGs)

    summary = np.clip(summary, a_min= min_clip, a_max= max_clip)

    # summary[np.abs(summary) < 1e-14] = 0

    return summary
