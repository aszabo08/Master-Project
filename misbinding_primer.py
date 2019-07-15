


from PCR_project import *





# amplicon length - primer length - n = 1000 - 15 -10 = 975         -------->  the longest extended length is 975 -1, when the product is only one nucleotide short
#                                                                              the shortest extended length is 0
# average extended length = (974 + 0)/2 = 487

misbinding_extended_length = 487

length_of_L = 487 + 15 + 10         # 512



Tm_misbinding_primer = 303.15               # 30 degrees

Tm_misbinding_extended_primer = 308.15      # 35 degree

Tm_misbinding_single_substrate = 343.15     # 70 degree    L  ( length of the created shorter product is 512 nt ---->  Neb calculator )

Tm_misbinding_double_substrate = 345.15     # 72 degree    L1L2  -----> it is higher as assumingly the there are no mismatches in the complex while in S1L2 there could be more mismatched nt pairs due to the misbinding



def misbinding_primer_1(values, t, T, dGs):



    S1 = values[1]
    P2 = values[4]
    S1M2 = values[17]


    kf6 = forward_rate

    #exponent_2a = dGs[1]/(R*T)

    #exponent_2a = np.clip((dGs[1]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_6a = exponent_clipping(dGs[4]/(R*T))

    kr6a = kf6 * np.exp(exponent_6a)

    kf6, kr6a = clipping(kf6, kr6a)



    #rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2

    #rate_S1P2_bind = np.clip((kf2 * S1 * P2 - kr2a * S1P2), a_min= min_clip, a_max= max_clip)

    rate_S1M2_bind = rate_clipping(kf6 * S1 * P2 - kr6a * S1M2)


    y = np.zeros(33)


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

    y = np.zeros(33)

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

    y = np.zeros(33)

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


    y = np.zeros(33)


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

    exponent_8 = exponent_clipping(dGs[6]/(R*T))

    kr8 = kf8 * np.exp(exponent_8)


    kf8, kr8 = clipping(kf8, kr8)



    #rate_den = - kr1 * S1S2 + kf1 * S1 * S2

    #rate_den = np.clip((- kr1 * S1S2 + kf1 * S1 * S2), a_min= min_clip, a_max= max_clip)

    rate_den = rate_clipping(- kr8 * S1L2 + kf8 * S1 * L2)


    y = np.zeros(33)

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

    y = np.zeros(33)

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
    Q2 = values[12]
    S1N2 = values[22]



    kf10 = forward_rate

    #exponent_5a = dGs[2]/(R*T)

    #exponent_5a = np.clip((dGs[2]/(R*T)), a_min= None, a_max= max_exponent)

    exponent_10a = exponent_clipping(dGs[5]/(R*T))

    kr10a = kf10 * np.exp(exponent_10a)

    kf10, kr10a = clipping(kf10, kr10a)

    #rate_S1Q2_bind = kf5 * S1 * Q2 - kr5a * S1Q2

    rate_S1N2_bind = rate_clipping(kf10 * S1 * Q2 - kr10a * S1N2)

    y = np.zeros(33)


    y[1] = - rate_S1N2_bind
    y[12] = - rate_S1N2_bind
    y[22] = rate_S1N2_bind



    return y















def short_misbinding_primer(values, t, T, dGs):


    L2P1 = values[23]
    P1 = values[3]
    L2 = values[21]


    L1P2 = values[29]
    P2 = values[4]
    L1 = values[28]



    kf11 = forward_rate


    #exponent_1 = dGs[0]/(R*T)

    exponent_11a = exponent_clipping(dGs[1]/(R*T))

    kr11a = kf11 * np.exp(exponent_11a)


    kf11, kr11a = clipping(kf11, kr11a)

    rate_L2P1 = rate_clipping(- kr11a * L2P1 + kf11 * P1 * L2)



    #exponent_1 = dGs[0]/(R*T)

    exponent_11b = exponent_clipping(dGs[0]/(R*T))

    kr11b = kf11 * np.exp(exponent_11b)


    kf11, kr11b = clipping(kf11, kr11b)

    rate_L1P2 = rate_clipping(- kr11b * L1P2 + kf11 * P2 * L1)



    y = np.zeros(33)

    y[23] = rate_L2P1        # concentration of P1L2

    y[3] = -rate_L2P1        # concentration of P1

    y[21] = -rate_L2P1       # concentration of L2





    y[29] = rate_L1P2        # concentration of L1P2

    y[4] = -rate_L1P2        # concentration of P2

    y[28] = -rate_L1P2      # concentration of L1



    return y




def short_misbinding_polymerase_1(values, t, T, dGs):


    L2P1 = values[23]
    E = values[7]
    L2P1E = values[24]



    L1P2 = values[29]
    #E = values[7]
    L1P2E = values[30]



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


    y = np.zeros(33)

    y[23] = - rate_poly_L2P1_bind
    #y[7] = - rate_poly_L2P1_bind
    y[24] = rate_poly_L2P1_bind


    y[29] = - rate_poly_L1P2_bind
    y[7] = enzyme_misbinding
    y[30] = rate_poly_L1P2_bind



    return y



def short_misbinding_primer_ext_1(values, t, T, dGs):

    """ The primer is extended by a few nucleotides to ensure the primer binding
        to the substrate without dissociation when reaching the extension temperature

        In this model 99.5 % of the primer - substrate complexes stay together, while 0.5 % of them will melt
        at the extension temperature"""


    L2P1E = values[24]
    dNTP = values[10]


    L1P2E = values[30]




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

    y = np.zeros(33)

    y[24] = - rate_ext_1                 # concentration of L2P1E
    #y[10] = - n * rate_ext_1             # concentration of dNTP
    y[25] = rate_ext_1                   # concentration of L2Q1E


    y[30] = - rate_ext_2                 # concentration of L1P2E
    y[10] = misbinding_nucleotide            # concentration of dNTP
    y[31] = rate_ext_2                   # concentration of L1Q2E


    return y




def short_misbinding_primer_ext_2(values, t, T, dGs):



    L2Q1E = values[25]
    dNTP = values[10]


    L1Q2E = values[31]



    ce_L = taq_nt_per_s(T)

    # reaction: S1N2E + misbinding_extended_length * dNTP ---> L1L2 + E

    #rate_ext_Q1 = (ce_Q / extended_length) * S1Q2E * dNTP

    #rate_ext_Q1 = np.clip(((ce_Q / extended_length) * S1Q2E * dNTP), a_min= min_clip, a_max= max_clip)

    rate_ext_QL_1 = rate_clipping((ce_L / misbinding_extended_length) * L2Q1E * dNTP)




    rate_ext_QL_2 = rate_clipping((ce_L / misbinding_extended_length) * L1Q2E * dNTP)


    mis_product = rate_clipping(rate_ext_QL_1 + rate_ext_QL_2)

    mis_nucleotide = rate_clipping(- misbinding_extended_length * rate_ext_QL_1 - misbinding_extended_length * rate_ext_QL_2)


    y = np.zeros(33)


    #y[7] = rate_ext_QL_1          # concentration of enzyme
    #y[10] = - misbinding_extended_length * rate_ext_QL_1           # concentration of dNTP
    y[25] = - rate_ext_QL_1       # concentration of L2Q1E
    #y[27] = rate_ext_QL_1         # concentration of L1L2



    y[7] = mis_product          # concentration of enzyme
    y[10] = mis_nucleotide          # concentration of dNTP
    y[31] = - rate_ext_QL_2      # concentration of L1Q2E
    y[26] = mis_product        # concentration of L1L2


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
    L2Q1E = values[25]
    L2Q1 = values[27]


    L1Q2E = values[31]
    L1Q2 = values[32]



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

    mis_enzyme = rate_clipping(- rate_poly_Q1L2_bind - rate_poly_L1Q2_bind)



    y = np.zeros(33)

    y[27] = - rate_poly_Q1L2_bind
    #y[7] = - rate_poly_Q1L2_bind
    y[25] = rate_poly_Q1L2_bind



    y[32] = - rate_poly_L1Q2_bind
    #y[7] = - rate_poly_L1Q2_bind
    y[31] = rate_poly_L1Q2_bind

    y[7] = mis_enzyme




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
    L2Q1 = values[27]



    Q2 = values[12]
    L1 = values[28]
    L1Q2 = values[32]



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

    y = np.zeros(33)


    y[11] = - rate_Q1L2_bind
    y[21] = - rate_Q1L2_bind
    y[27] = rate_Q1L2_bind





    y[12] = - rate_L1Q2_bind
    y[28] = - rate_L1Q2_bind
    y[32] = rate_L1Q2_bind





    return y






def L_misbinding_denaturation(values, t, T, dGs):


    L1L2 = values[26]
    L1 = values[28]
    L2 = values[21]


    kf15 = forward_rate


    #exponent_1 = dGs[0]/(R*T)

    exponent_15 = exponent_clipping(dGs[7]/(R*T))

    kr15 = kf15 * np.exp(exponent_15)


    kf15, kr15 = clipping(kf15, kr15)

    rate_den = rate_clipping(- kr15 * L1L2 + kf15 * L1 * L2)


    y = np.zeros(33)

    y[26] = rate_den        # concentration of L1L2

    y[28] = -rate_den        # concentration of L1

    y[21] = -rate_den       # concentration of L2



    return y







def PCR_reaction_with_misbinding(values, t, T, dGs):  # not using rate clipping cos the array data

    summary = denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + primer_binding_2(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + primer_ext_2(values, t, T, dGs) + taq_denaturation(values, t, T, dGs) + misbinding_primer_1(values, t, T, dGs) + misbinding_polymerase_1(values, t, T, dGs) + misbinding_primer_ext_1(values, t, T, dGs) + misbinding_primer_ext_2(values, t, T, dGs) + misbinding_denaturation(values, t, T, dGs) + misbinding_polymerase_2(values, t, T, dGs) + misbinding_primer_2(values, t, T, dGs) + short_misbinding_primer(values, t, T, dGs) + short_misbinding_polymerase_1(values, t, T, dGs) + short_misbinding_primer_ext_1(values, t, T, dGs) + short_misbinding_primer_ext_2(values, t, T, dGs) + short_misbinding_polymerase_2(values, t, T, dGs) + short_misbinding_primer_2(values, t, T, dGs) + L_misbinding_denaturation(values, t, T, dGs)

    summary = np.clip(summary, a_min= min_clip, a_max= max_clip)

    # summary[np.abs(summary) < 1e-14] = 0

    return summary


























def misbinding_only_one_integration(values, number):

    number = number - 1


    functions_with_misbinding = [denaturation, primer_binding_1, polymerase_binding_1, primer_ext_1, polymerase_binding_2, primer_binding_2, primer_ext_2, PCR_reaction, misbinding_primer_1, misbinding_polymerase_1, misbinding_primer_ext_1, misbinding_primer_ext_2, misbinding_denaturation, misbinding_polymerase_2, misbinding_primer_2, short_misbinding_primer, short_misbinding_polymerase_1, short_misbinding_primer_ext_1, short_misbinding_primer_ext_2, short_misbinding_polymerase_2, short_misbinding_primer_2, L_misbinding_denaturation, PCR_reaction_with_misbinding]

    functions_name_with_misbinding = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer_ext_2", "PCR_reaction",  "misbinding_primer_1", "misbinding_polymerase_1", "misbinding_primer_ext_1", "misbinding_primer_ext_2", "misbinding_denaturation", "misbinding_polymerase_2", "misbinding_primer_2", "short_misbinding_primer", "short_misbinding_polymerase_1", "short_misbinding_primer_ext_1", "short_misbinding_primer_ext_2", "short_misbinding_polymerase_2", "short_misbinding_primer_2", "L_misbinding_denaturation", "PCR_reaction_with_misbinding"]



    #no PCR, no PCR misbinding
    #all_function_names = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer_ext_2", "PCR_reaction",  "misbinding_primer_1", "misbinding_polymerase_1", "misbinding_primer_ext_1", "misbinding_primer_ext_2", "misbinding_denaturation", "misbinding_polymerase_2", "misbinding_primer_2", "short_misbinding_primer", "short_misbinding_polymerase_1", "short_misbinding_primer_ext_1", "short_misbinding_primer_ext_2", "short_misbinding_polymerase_2", "short_misbinding_primer_2", "L_misbinding_denaturation", "PCR_reaction_with_misbinding"]




    concentration = np.empty((number_time_points, 33))

    dGs = [0 for x in range(8)]



    #before_nt_number = all_nucleotide(values)

    #all_umol_1 = u_molar_concentration(values)


    for i in range(number_cycles):


        dGs[0] = (Tm_S1S2 - Tden) * dS

        dGs[1] = (Tm_primer - Tden) * dS

        dGs[2] = (Tm_extended_primer - Tden) * dS

        dGs[3] = (Tm_enzyme - Tden) * dS

        dGs[4] = (Tm_misbinding_primer - Tden) * dS

        dGs[5] = (Tm_misbinding_extended_primer - Tden) * dS

        dGs[6] = (Tm_misbinding_single_substrate - Tden) * dS

        dGs[7] = (Tm_misbinding_double_substrate - Tden) * dS


        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        #print("dGs_den", dGs)




        integration_den = odeint(functions_with_misbinding[number], values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs))

        concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den

        dGs[0] = (Tm_S1S2 - Tanneal) * dS

        dGs[1] = (Tm_primer - Tanneal) * dS

        dGs[2] = (Tm_extended_primer - Tanneal) * dS

        dGs[3] = (Tm_enzyme - Tanneal) * dS

        dGs[4] = (Tm_misbinding_primer - Tanneal) * dS

        dGs[5] = (Tm_misbinding_extended_primer - Tanneal) * dS

        dGs[6] = (Tm_misbinding_single_substrate - Tanneal) * dS

        dGs[7] = (Tm_misbinding_double_substrate - Tanneal) * dS



        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        #print("dGs_anneals", dGs)

        integration_anneal = odeint(functions_with_misbinding[number], integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs))

        concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal

        dGs[0] = (Tm_S1S2 - Text) * dS

        dGs[1] = (Tm_primer - Text) * dS

        dGs[2] = (Tm_extended_primer - Text) * dS

        dGs[3] = (Tm_enzyme - Text) * dS

        dGs[4] = (Tm_misbinding_primer - Text) * dS

        dGs[5] = (Tm_misbinding_extended_primer - Text) * dS

        dGs[6] = (Tm_misbinding_single_substrate - Text) * dS

        dGs[7] = (Tm_misbinding_double_substrate - Text) * dS

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        #print("dGs_text", dGs)

        integration_ext = odeint(functions_with_misbinding[number], integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs))

        concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


        values = integration_ext[-1]




    print("The concentration of the 17 species at the end of", functions_name_with_misbinding[number], ":", values)

    #after_nt_number = all_nucleotide(values)

    #all_umol_2 = u_molar_concentration(values)

    #print("The difference in nt number after", functions_plus_name[number], ":", before_nt_number - after_nt_number, "\n")

    #print("The difference in micromolar concentration after", functions_plus_name[number], ":", all_umol_1 - all_umol_2, "\n")


    #Plotting the concentrations over time

    plt.figure(1)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)

    #plots1 = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]

    plots1 = []

    for i in range(10):

        plots1.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(10):

        plt.subplot(2, 5, i+1)

        plt.gca().set_title(new_species[plots1[i]])


        plt.plot(time, concentration[:, plots1[i]])

        #plt.ylim([0, y_top_limit[i]])


        #plt.legend([new_species[plots1[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time (s)")
        plt.ylabel("Concentration (uM)")
        plt.tight_layout(w_pad=0.8, h_pad=0.8)


    plt.figure(2)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)



    plots2 = []

    for i in range(10, 17):

        plots2.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(7):

        plt.subplot(2, 4, i+1)

        plt.gca().set_title(new_species[plots2[i]])

        plt.plot(time, concentration[:, plots2[i]])

        #plt.ylim([0, y_top_limit[i]])


        #plt.legend([new_species[plots2[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time (s)")
        plt.ylabel("Concentration (uM)")
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.35, hspace=0.4)



    plt.figure(3)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)

    plots3 = []

    for i in range(17, 27):

        plots3.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(10):

        plt.subplot(2, 5, i+1)

        plt.gca().set_title(new_species[plots3[i]])

        plt.plot(time, concentration[:, plots3[i]])

        #plt.ylim([0, y_top_limit[i]])


        #plt.legend([new_species[plots3[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time (s)")
        plt.ylabel("Concentration (uM)")
        plt.tight_layout()



    plt.figure(4)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)


    plots4 = []

    for i in range(27, 33):

        plots4.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(6):

        plt.subplot(2, 4, i+1)

        plt.gca().set_title(new_species[plots4[i]])

        plt.plot(time, concentration[:, plots4[i]])

        #plt.ylim([0, y_top_limit[i]])


        #plt.legend([new_species[plots4[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time (s)")
        plt.ylabel("Concentration (uM)")
        plt.tight_layout()






    plt.show()




    return values


if __name__ == '__main__':

    print(values)

    misbinding_only_one_integration(values, 23)    # pcr with misbinding
