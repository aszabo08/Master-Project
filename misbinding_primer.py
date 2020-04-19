


from PCR_project import *
from data_input import *



##########  Lengths and melting temperatures of misbinding species  ##########


# amplicon length - primer length - n = 1000 - 15 -10 = 975         -------->  the longest extended length is 975 -1, when the product is only one nucleotide short
#                                                                              the shortest extended length is 0
# average extended length = (974 + 0)/2 = 487

misbinding_extended_length = 487


# for the experiment
# misbinding_extended_length = 1086


#length_of_L = 487 + 15 + 10         # 512

length_of_L = misbinding_extended_length + primer_length + n

#mismatch = int(primer_length * 0.75)

mismatch = int(0.5 * primer_length)

#mismatch = int(0.8 * primer_length)

length_misbinding_primer = primer_length - mismatch

length_misbinding_extended_primer = extended_primer - mismatch

length_misbinding_single_substrate = length_of_L - mismatch



Tm_misbinding_primer = (Tmax * length_misbinding_primer * dH) / ((length_misbinding_primer + K) * dS)                                       # 342.71 K

Tm_misbinding_extended_primer = (Tmax * length_misbinding_extended_primer * dH) / ((length_misbinding_extended_primer + K) * dS)            # 346.87 K

Tm_misbinding_single_substrate = (Tmax * length_misbinding_single_substrate * dH) / ((length_misbinding_single_substrate + K) * dS)         # 356.168 K

Tm_misbinding_double_substrate = (Tmax * length_of_L * dH) / ((length_of_L + K) * dS)                                                       # 356.17 K



##########  Interactions of PCR with primer misbinding  ##########



def misbinding_primer_1(values, t, T, dGs):

    S1 = values[1]
    P2 = values[4]
    S1M2 = values[17]


    kf6 = forward_rate

    exponent_6a = exponent_clipping(dGs[4]/(R*T))

    kr6a = kf6 * np.exp(exponent_6a)

    kf6, kr6a = clipping(kf6, kr6a)

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

    exponent_7a = exponent_clipping(dGs[3]/(R*T))

    kr7a = kf7 * np.exp(exponent_7a)

    kf7, kr7a = clipping(kf7, kr7a)

    rate_poly_S1M2_bind = rate_clipping(kf7 * S1M2 * E - kr7a * S1M2E)

    y = np.zeros(33)

    y[17] = - rate_poly_S1M2_bind
    y[7] = - rate_poly_S1M2_bind
    y[18] = rate_poly_S1M2_bind



    return y





def misbinding_primer_ext_1(values, t, T, dGs):



    S1M2E = values[18]
    dNTP = values[10]


    ce = polymerase_nt_per_s(T)


    rate_ext_1 = rate_clipping((ce / n) * S1M2E * dNTP)

    y = np.zeros(33)

    y[18] = - rate_ext_1  # concentration of S1M2E
    y[10] = - n * rate_ext_1  # concentration of dNTP
    y[19] = rate_ext_1  # concentration of S1N2E


    return y




def misbinding_primer_ext_2(values, t, T, dGs):



    S1N2E = values[19]
    dNTP = values[10]

    ce_N = polymerase_nt_per_s(T)

    # reaction: S1N2E + misbinding_extended_length * dNTP ---> L1L2 + E

    rate_ext_N1 = rate_clipping((ce_N / misbinding_extended_length) * S1N2E * dNTP)


    y = np.zeros(33)


    y[7] = rate_ext_N1
    y[10] = - misbinding_extended_length * rate_ext_N1           # concentration of dNTP
    y[19] = - rate_ext_N1       # concentration of S1N2E
    y[20] = rate_ext_N1         # concentration of S1L2

    return y




def misbinding_denaturation(values, t, T, dGs):


    S1L2 = values[20]
    S1 = values[1]
    L2 = values[21]

    kf8 = forward_rate

    exponent_8 = exponent_clipping(dGs[6]/(R*T))

    kr8 = kf8 * np.exp(exponent_8)


    kf8, kr8 = clipping(kf8, kr8)

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

    exponent_9a = exponent_clipping(dGs[3]/(R*T))


    kr9a = kf9 * np.exp(exponent_9a)

    kf9, kr9a = clipping(kf9, kr9a)

    rate_poly_S1N2_bind = rate_clipping(kf9 * S1N2 * E - kr9a * S1N2E)

    y = np.zeros(33)

    y[22] = - rate_poly_S1N2_bind
    y[7] = - rate_poly_S1N2_bind
    y[19] = rate_poly_S1N2_bind



    return y





def misbinding_primer_2(values, t, T, dGs):




    S1 = values[1]
    Q2 = values[12]
    S1N2 = values[22]



    kf10 = forward_rate

    exponent_10a = exponent_clipping(dGs[5]/(R*T))

    kr10a = kf10 * np.exp(exponent_10a)

    kf10, kr10a = clipping(kf10, kr10a)


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


    exponent_11a = exponent_clipping(dGs[1]/(R*T))

    kr11a = kf11 * np.exp(exponent_11a)


    kf11, kr11a = clipping(kf11, kr11a)

    rate_L2P1 = rate_clipping(- kr11a * L2P1 + kf11 * P1 * L2)


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

    exponent_12a = exponent_clipping(dGs[3]/(R*T))

    kr12a = kf12 * np.exp(exponent_12a)

    kf12, kr12a = clipping(kf12, kr12a)

    rate_poly_L2P1_bind = rate_clipping(kf12 * L2P1 * E - kr12a * L2P1E)


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


    L2P1E = values[24]
    dNTP = values[10]


    L1P2E = values[30]

    ce = polymerase_nt_per_s(T)

    rate_ext_1 = rate_clipping((ce / n) * L2P1E * dNTP)


    rate_ext_2 = rate_clipping((ce / n) * L1P2E * dNTP)


    misbinding_nucleotide = rate_clipping(- n * rate_ext_1 - n * rate_ext_2)

    y = np.zeros(33)

    y[24] = - rate_ext_1                 # concentration of L2P1E
    y[25] = rate_ext_1                   # concentration of L2Q1E


    y[30] = - rate_ext_2                 # concentration of L1P2E
    y[10] = misbinding_nucleotide            # concentration of dNTP
    y[31] = rate_ext_2                   # concentration of L1Q2E


    return y




def short_misbinding_primer_ext_2(values, t, T, dGs):



    L2Q1E = values[25]
    dNTP = values[10]


    L1Q2E = values[31]



    ce_L = polymerase_nt_per_s(T)

    # reaction: S1N2E + misbinding_extended_length * dNTP ---> L1L2 + E

    rate_ext_QL_1 = rate_clipping((ce_L / misbinding_extended_length) * L2Q1E * dNTP)

    rate_ext_QL_2 = rate_clipping((ce_L / misbinding_extended_length) * L1Q2E * dNTP)


    mis_product = rate_clipping(rate_ext_QL_1 + rate_ext_QL_2)

    mis_nucleotide = rate_clipping(- misbinding_extended_length * rate_ext_QL_1 - misbinding_extended_length * rate_ext_QL_2)


    y = np.zeros(33)


    y[25] = - rate_ext_QL_1       # concentration of L2Q1E



    y[7] = mis_product          # concentration of enzyme
    y[10] = mis_nucleotide          # concentration of dNTP
    y[31] = - rate_ext_QL_2      # concentration of L1Q2E
    y[26] = mis_product        # concentration of L1L2


    return y



def short_misbinding_polymerase_2(values, t, T, dGs):



    E = values[7]
    L2Q1E = values[25]
    L2Q1 = values[27]


    L1Q2E = values[31]
    L1Q2 = values[32]



    kf13 = forward_rate

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
    y[25] = rate_poly_Q1L2_bind



    y[32] = - rate_poly_L1Q2_bind
    y[31] = rate_poly_L1Q2_bind

    y[7] = mis_enzyme




    return y




def short_misbinding_primer_2(values, t, T, dGs):


    Q1 = values[11]
    L2 = values[21]
    L2Q1 = values[27]

    Q2 = values[12]
    L1 = values[28]
    L1Q2 = values[32]

    kf14 = forward_rate

    exponent_14a = exponent_clipping(dGs[2]/(R*T))

    kr14a = kf14 * np.exp(exponent_14a)

    kf14, kr14a = clipping(kf14, kr14a)

    rate_Q1L2_bind = rate_clipping(kf14 * Q1 * L2 - kr14a * L2Q1)

    exponent_14b = exponent_clipping(dGs[2]/(R*T))

    kr14b = kf14 * np.exp(exponent_14b)

    kf14, kr14b = clipping(kf14, kr14b)

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

    exponent_15 = exponent_clipping(dGs[7]/(R*T))

    kr15 = kf15 * np.exp(exponent_15)


    kf15, kr15 = clipping(kf15, kr15)

    rate_den = rate_clipping(- kr15 * L1L2 + kf15 * L1 * L2)


    y = np.zeros(33)

    y[26] = rate_den        # concentration of L1L2

    y[28] = -rate_den        # concentration of L1

    y[21] = -rate_den       # concentration of L2



    return y







def PCR_reaction_with_misbinding(values, t, T, dGs):

    """ This functions adds up all the interactions included in PCR with primer misbinding """

    summary = denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs) + primer_binding_2(values, t, T, dGs) + polymerase_binding_1(values, t, T, dGs) + polymerase_binding_2(values, t, T, dGs) + primer_ext_1(values, t, T, dGs) + primer_ext_2(values, t, T, dGs) + enzyme_denaturation(values, t, T, dGs) + misbinding_primer_1(values, t, T, dGs) + misbinding_polymerase_1(values, t, T, dGs) + misbinding_primer_ext_1(values, t, T, dGs) + misbinding_primer_ext_2(values, t, T, dGs) + misbinding_denaturation(values, t, T, dGs) + misbinding_polymerase_2(values, t, T, dGs) + misbinding_primer_2(values, t, T, dGs) + short_misbinding_primer(values, t, T, dGs) + short_misbinding_polymerase_1(values, t, T, dGs) + short_misbinding_primer_ext_1(values, t, T, dGs) + short_misbinding_primer_ext_2(values, t, T, dGs) + short_misbinding_polymerase_2(values, t, T, dGs) + short_misbinding_primer_2(values, t, T, dGs) + L_misbinding_denaturation(values, t, T, dGs)

    summary = np.clip(summary, a_min= min_clip, a_max= max_clip)

    return summary



##########  PCR with misbinding integrations and assessments ##########



def PCR_misbinding_integration(values):

    """ This function integrates the PCR reaction with primer misbinding over t time """



    concentration = np.empty((number_time_points, 33))

    dGs = [0 for x in range(8)]

    dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_initial_den * dS)

    dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_initial_den * dS)

    dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_initial_den * dS)

    dGs[3] = (Tm_enzyme - T_initial_den) * dS

    dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (T_initial_den * dS)

    dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (T_initial_den * dS)

    dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (T_initial_den * dS)

    dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (T_initial_den * dS)

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

        dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Tden * dS)

        dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Tden * dS)

        dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Tden * dS)

        dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Tden * dS)

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        # integration_den = odeint(PCR_reaction_with_misbinding, values, time[(t_initial_den * steps -1 + total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)], args=(Tden, dGs),  mxstep=5000000)
        #
        # concentration[t_initial_den * steps -1 + (total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)] = integration_den


        integration_den = odeint(PCR_reaction_with_misbinding, values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

        concentration[(total * i * steps): ((total * i + tden) * steps)] = integration_den



        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tanneal * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tanneal * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tanneal * dS)

        dGs[3] = (Tm_enzyme - Tanneal) * dS

        dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Tanneal * dS)

        dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Tanneal * dS)

        dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Tanneal * dS)

        dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Tanneal * dS)


        dGs = np.clip(dGs, a_min=None, a_max=1e+12)


        # integration_anneal = odeint(PCR_reaction_with_misbinding, integration_den[-1], time[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs),  mxstep=5000000)

        # concentration[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)] = integration_anneal


        integration_anneal = odeint(PCR_reaction_with_misbinding, integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)


        concentration[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal



        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Text * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Text * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Text * dS)

        dGs[3] = (Tm_enzyme - Text) * dS

        dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Text * dS)

        dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Text * dS)

        dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Text * dS)

        dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Text * dS)

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        # integration_ext = odeint(PCR_reaction_with_misbinding, integration_anneal[-1], time[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)], args=(Text, dGs),  mxstep=5000000)
        #
        # concentration[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)] = integration_ext


        integration_ext = odeint(PCR_reaction_with_misbinding, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)

        concentration[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext


        values = integration_ext[-1]


    dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_cooling_down * dS)

    dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_cooling_down  * dS)

    dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_cooling_down  * dS)

    dGs[3] = (Tm_enzyme - T_cooling_down ) * dS

    dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (T_cooling_down  * dS)

    dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (T_cooling_down  * dS)

    dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (T_cooling_down * dS)

    dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (T_cooling_down  * dS)

    dGs = np.clip(dGs, a_min=None, a_max=1e+12)

    # integration_cool = odeint(PCR_reaction_with_misbinding, integration_ext[-1], time[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)], args=(T_cooling_down, dGs), mxstep=5000000)
    #
    # concentration[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)] = integration_cool[-1]

    #values = integration_cool[-1]

    print("The concentration of the species at the end of PCR misbinding integration:", values)

    print("The concentration of S1S2 in ng/ul:", uM_to_ng_per_ul(values[0], amplicon_length))
    print("The concentration of L1L2 in ng/ul:", uM_to_ng_per_ul(values[26], length_of_L))

    #print("total misbinding conc", misbinding_PCR_total_concentration(concentration, time))

    purity_over_total_yield(concentration, time)


    #Plotting the concentrations over time

    plt.figure(1)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 18,  fontweight = 'bold', y = 0.95)

    #plots1 = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]

    plots1 = []

    for i in range(10):

        plots1.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(10):

        plt.subplot(2, 5, i+1)

        plt.gca().set_title(new_species[plots1[i]], fontweight = 'bold')


        plt.plot(time, concentration[:, plots1[i]])

        #plt.ylim([0, y_top_limit[i]])

        plt.xlabel("Time (s)", fontsize = 12)
        plt.ylabel("Concentration (uM)", fontsize = 12)
        plt.subplots_adjust(wspace = 0.38)


    plt.figure(2)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 18,  fontweight = 'bold', y = 0.95)

    plots2 = []

    for i in range(10, 17):

        plots2.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(7):

        plt.subplot(2, 4, i+1)

        plt.gca().set_title(new_species[plots2[i]], fontweight = 'bold')

        plt.plot(time, concentration[:, plots2[i]])

        #plt.ylim([0, y_top_limit[i]])

        plt.xlabel("Time (s)", fontsize = 12)
        plt.ylabel("Concentration (uM)", fontsize = 12)
        plt.subplots_adjust(wspace = 0.38)



    plt.figure(3)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 18,  fontweight = 'bold', y = 0.95)

    plots3 = []

    for i in range(17, 27):

        plots3.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(10):

        plt.subplot(2, 5, i+1)

        plt.gca().set_title(new_species[plots3[i]], fontweight = 'bold')

        plt.plot(time, concentration[:, plots3[i]])

        #plt.ylim([0, y_top_limit[i]])

        plt.xlabel("Time (s)", fontsize = 12)
        plt.ylabel("Concentration (uM)", fontsize = 12)
        plt.subplots_adjust(wspace = 0.38)



    plt.figure(4)

    plt.suptitle("Change of the species' concentrations over time" , fontsize = 18,  fontweight = 'bold', y = 0.95)


    plots4 = []

    for i in range(27, 33):

        plots4.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(6):

        plt.subplot(2, 4, i+1)

        plt.gca().set_title(new_species[plots4[i]], fontweight = 'bold')

        plt.plot(time, concentration[:, plots4[i]])

        #plt.ylim([0, y_top_limit[i]])


        plt.xlabel("Time (s)", fontsize = 12)
        plt.ylabel("Concentration (uM)", fontsize = 12)
        plt.subplots_adjust(wspace = 0.38)

    plt.show()

    return values





def extended_PCR_misbinding_integration(values):

    """ This function integrates the PCR reaction with primer misbinding with added initial denaturation and final elongation"""


    concentration = np.empty((number_time_points, 33))

    dGs = [0 for x in range(8)]


    dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_initial_den * dS)

    dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_initial_den * dS)

    dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_initial_den * dS)

    dGs[3] = (Tm_enzyme - T_initial_den) * dS

    dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (T_initial_den * dS)

    dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (T_initial_den * dS)

    dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (T_initial_den * dS)

    dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (T_initial_den * dS)

    dGs = np.clip(dGs, a_min=None, a_max=1e+12)


    integration_initial_den = odeint(PCR_reaction, values, time[0: t_initial_den * steps], args=(T_initial_den, dGs), mxstep=5000000)

    concentration[0: t_initial_den * steps] = integration_initial_den

    values = integration_initial_den[-1]


    for i in range(number_cycles):

        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tden * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tden * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tden * dS)

        dGs[3] = (Tm_enzyme - Tden) * dS

        dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Tden * dS)

        dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Tden * dS)

        dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Tden * dS)

        dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Tden * dS)

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        integration_den = odeint(PCR_reaction_with_misbinding, values, time[(t_initial_den * steps -1 + total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)], args=(Tden, dGs),  mxstep=5000000)

        concentration[t_initial_den * steps -1 + (total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)] = integration_den


        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tanneal * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tanneal * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tanneal * dS)

        dGs[3] = (Tm_enzyme - Tanneal) * dS

        dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Tanneal * dS)

        dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Tanneal * dS)

        dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Tanneal * dS)

        dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Tanneal * dS)


        dGs = np.clip(dGs, a_min=None, a_max=1e+12)


        integration_anneal = odeint(PCR_reaction_with_misbinding, integration_den[-1], time[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs),  mxstep=5000000)

        concentration[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)] = integration_anneal


        dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Text * dS)

        dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Text * dS)

        dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Text * dS)

        dGs[3] = (Tm_enzyme - Text) * dS

        dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Text * dS)

        dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Text * dS)

        dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Text * dS)

        dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Text * dS)

        dGs = np.clip(dGs, a_min=None, a_max=1e+12)

        integration_ext = odeint(PCR_reaction_with_misbinding, integration_anneal[-1], time[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)], args=(Text, dGs),  mxstep=5000000)

        concentration[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)] = integration_ext


        values = integration_ext[-1]


    dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_cooling_down * dS)

    dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_cooling_down  * dS)

    dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_cooling_down  * dS)

    dGs[3] = (Tm_enzyme - T_cooling_down ) * dS

    dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (T_cooling_down  * dS)

    dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (T_cooling_down  * dS)

    dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (T_cooling_down * dS)

    dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (T_cooling_down  * dS)

    dGs = np.clip(dGs, a_min=None, a_max=1e+12)


    integration_cool = odeint(PCR_reaction_with_misbinding, integration_ext[-1], time[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)], args=(T_cooling_down, dGs), mxstep=5000000)

    concentration[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)] = integration_cool[-1]

    values = integration_cool[-1]

    print("The concentration of the 17 species at the end of PCR misbinding integration:", values)

    print("The concentration of S1S2 in ng/ul:", uM_to_ng_per_ul(values[0], amplicon_length))
    print("The concentration of L1L2 in ng/ul:", uM_to_ng_per_ul(values[26], length_of_L))

    #print("total misbinding conc", misbinding_PCR_total_concentration(concentration, time))

    purity_over_total_yield(concentration, time)


    #Plotting the concentrations over time
    plt.figure(1)

    plt.suptitle("Change of the species' concentrations over time", fontweight = 'bold')

    #plots1 = [0, 1, 3, 5, 7, 8, 10, 11, 13, 15]

    plots1 = []

    for i in range(10):

        plots1.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(10):

        plt.subplot(2, 5, i+1)

        plt.gca().set_title(new_species[plots1[i]], fontweight = 'bold')
        plt.gca().title.set_position([0.5, 1.05])


        plt.plot(time, concentration[:, plots1[i]])

        #plt.ylim([0, y_top_limit[i]])

        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,2))
        plt.xlabel("Time (s)")
        plt.ylabel("Concentration (uM)")


    plt.subplots_adjust(hspace=0.5, wspace=1)


    plt.figure(2)

    plt.suptitle("Change of the species' concentrations over time", fontweight = 'bold')



    plots2 = []

    for i in range(10, 17):

        plots2.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(7):

        plt.subplot(2, 4, i+1)

        plt.gca().set_title(new_species[plots2[i]], fontweight = 'bold')
        plt.gca().title.set_position([0.5, 1.05])

        plt.plot(time, concentration[:, plots2[i]])

        #plt.ylim([0, y_top_limit[i]])

        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,2))
        plt.xlabel("Time (s)")
        plt.ylabel("Concentration (uM)")



    plt.subplots_adjust(hspace=0.5, wspace=1)

    plt.figure(3)

    plt.suptitle("Change of the species' concentrations over time", fontweight = 'bold')

    plots3 = []

    for i in range(17, 27):

        plots3.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(10):

        plt.subplot(2, 5, i+1)

        plt.gca().set_title(new_species[plots3[i]], fontweight = 'bold')
        plt.gca().title.set_position([0.5, 1.05])

        plt.plot(time, concentration[:, plots3[i]])

        #plt.ylim([0, y_top_limit[i]])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,2))
        plt.xlabel("Time (s)")
        plt.ylabel("Concentration (uM)")

    plt.subplots_adjust(hspace=0.5, wspace=1)



    plt.figure(4)

    plt.suptitle("Change of the species' concentrations over time", fontweight = 'bold')


    plots4 = []

    for i in range(27, 33):

        plots4.append(i)

    y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(6):

        plt.subplot(2, 4, i+1)

        plt.gca().set_title(new_species[plots4[i]], fontweight = 'bold')
        plt.gca().title.set_position([0.5, 1.05])

        plt.plot(time, concentration[:, plots4[i]])

        #plt.ylim([0, y_top_limit[i]])
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,2))
        plt.xlabel("Time (s)")
        plt.ylabel("Concentration (uM)")


    plt.subplots_adjust(hspace=0.5, wspace=1)

    plt.show()

    return values




def misbinding_PCR_total_concentration(all_concentration, time_vector):

    """ This function calculates the total concentrations of the misbinding species throughout the PCR reaction with primer misbinding"""


    misbinding_single_species_PCR = ["S1", "S2", "P1", "P2", "Q1", "Q2", "E", "M2", "N2", "L1", "L2"]


    concentration_PCR = np.zeros((all_concentration.shape[0], len(misbinding_single_species_PCR)))         # the matrix s length is all time points

    indexes_of_species = [[] for i in range(len(misbinding_single_species_PCR))]


    for x in range(len(misbinding_single_species_PCR)):


        for i in new_species:

            if misbinding_single_species_PCR[x] in i:


                indexes_of_species[x].append(new_species.index(i))


    print(indexes_of_species)




    for x in range(all_concentration.shape[0]):


        for i in range(len(indexes_of_species)):

            summary_concentration = 0


            for n in range(len(indexes_of_species[i])):


                summary_concentration = summary_concentration + all_concentration[x, indexes_of_species[i][n]]



            concentration_PCR[x, i] = summary_concentration


    size = 18

    title_size = 22


    plt.figure(1)

    plt.suptitle("Total concentrations of single species over time with PCR misbinding" , fontsize = title_size)


    #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(len(misbinding_single_species_PCR)):

        plt.subplot(2, 6, i+1)


        plt.plot(time_vector, concentration_PCR[:, i])      #label = misbinding_single_species_PCR[i]

        plt.gca().set_title(misbinding_single_species_PCR[i])

        #plt.ylim([0, y_top_limit[i]])


        #plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})

        plt.xlabel("Time")
        plt.ylabel("Total concentration")




    plt.figure(2)

    plt.suptitle("The change of total concentrations of 4 single species\nwith primer misbinding" , fontsize=title_size,  fontweight= 'bold')

    highlighted_species = [0, 3, 5, 6]            # "S1", "S2", "P1", "P2", "Q1", "Q2", "E", "M2", "N2", "L1", "L2"


    #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(len(highlighted_species)):

        plt.subplot(2, 2, i+1)

        plt.tick_params(labelsize = 16)

        plt.gca().set_title(misbinding_single_species_PCR[highlighted_species[i]], fontweight = 'bold', fontsize = size)


        plt.plot(time_vector, concentration_PCR[:, highlighted_species[i]])   # label = misbinding_single_species_PCR[highlighted_species[i]]

        #plt.ylim([0, y_top_limit[i]])

        plt.xlabel("Time (s) ",  FontSize= size)
        plt.ylabel("Total concentration (uM)",  FontSize= size)


    plt.subplots_adjust(hspace = 0.3)



    plt.figure(3)

    plt.suptitle("The change of total concentrations of misbinding single species" , fontsize = title_size, fontweight= 'bold')

    highlighted_species = [7, 8, 10]            # "S1", "S2", "P1", "P2", "Q1", "Q2", "E", "M2", "N2", "L1", "L2"


    #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    for i in range(len(highlighted_species)):

        plt.subplot(2, 2, i+1)

        plt.tick_params(labelsize = 16)

        plt.gca().set_title(misbinding_single_species_PCR[highlighted_species[i]], fontweight = 'bold', fontsize = size )


        plt.plot(time_vector, concentration_PCR[:, highlighted_species[i]])  #label = misbinding_single_species_PCR[highlighted_species[i]]

        #plt.ylim([0, y_top_limit[i]])

        plt.xlabel("Time (s) ",  FontSize= size)
        plt.ylabel("Total concentration (uM)",  FontSize= size)

    plt.subplots_adjust(hspace = 0.3)

    plt.show()



    return concentration_PCR








first_concentration = [0 for i in range(33)]

#first_concentration[0] = 0.00015151515151515152

first_concentration[0] = 0.001

first_concentration[3], first_concentration[4] = 0.5, 0.5

first_concentration[7] = 0.2

first_concentration[10] = 200




second_concentration = [0 for i in range(33)]

#second_concentration[0] = 0.00015151515151515152

second_concentration[0] = 0.0001

second_concentration[3], second_concentration[4] = 0.5, 0.5

second_concentration[7] = 0.2

second_concentration[10] = 200


third_concentration = [0 for i in range(33)]

#third_concentration[0] = 0.00015151515151515152

third_concentration[0] = 0.00001

third_concentration[3], third_concentration[4] = 0.5, 0.5

third_concentration[7] = 0.2

third_concentration[10] = 200


fourth_concentration = [0 for i in range(33)]

fourth_concentration[0] = 0.00015151515151515152

fourth_concentration[3], third_concentration[4] = 0.5, 0.5

fourth_concentration[7] = 0.2

fourth_concentration[10] = 200


fifth_concentration = [0 for i in range(33)]

fifth_concentration[0] = 0.00015151515151515152

fifth_concentration[3], third_concentration[4] = 0.5, 0.5

fifth_concentration[7] = 0.2

fifth_concentration[10] = 200


#overall_concentration = [first_concentration, second_concentration, third_concentration, fourth_concentration, fifth_concentration]

overall_concentration = [first_concentration, second_concentration, third_concentration]

initial_overall = overall_concentration


T_multiple_annealing = [333.15, 334.15, 335.15, 336.15, 337.15]

multiple_mismatch = [int(0.25 * primer_length), int(0.5 * primer_length), int(0.75 * primer_length), int(0.85 * primer_length), int(0.90 * primer_length)]


label_s1s2 = ["S1S2 = " + str(first_concentration[0]), "S1S2 = " + str(second_concentration[0]), "S1S2 = " + str(third_concentration[0])]

#label_Text = ["Tanneal =" + str(T_multiple_annealing[0]), "Tanneal =" + str(T_multiple_annealing[1]), "Tanneal =" + str(T_multiple_annealing[2])]

label_mismatch = ["Length (nt) = " + str(multiple_mismatch[0]), "Length (nt) = " + str(multiple_mismatch[1]), "Length (nt) = " + str(multiple_mismatch[2]), "Length (nt) = " + str(multiple_mismatch[3]), "Length (nt) = " + str(multiple_mismatch[4])]





def purity_multiple_initial_conditions(overall_concentration):

    """ This function calculates the purity level of the final dsDNA
    under different settings """


    dGs = [0 for x in range(8)]


    concentration = np.zeros((len(overall_concentration), number_time_points, len(new_species)))


    for e in range(len(overall_concentration)):


        #Tanneal = T_multiple_annealing[e]

        # mismatch = multiple_mismatch[e]

        # length_misbinding_primer = primer_length - mismatch

        # length_misbinding_extended_primer = extended_primer - mismatch

        # length_misbinding_single_substrate = length_of_L - mismatch

        # Tm_misbinding_primer = (Tmax * length_misbinding_primer * dH) / ((length_misbinding_primer + K) * dS)                                       # 342.71 K

        # Tm_misbinding_extended_primer = (Tmax * length_misbinding_extended_primer * dH) / ((length_misbinding_extended_primer + K) * dS)            # 346.87 K

        # Tm_misbinding_single_substrate = (Tmax * length_misbinding_single_substrate * dH) / ((length_misbinding_single_substrate + K) * dS)         # 356.168 K

        # Tm_misbinding_double_substrate = (Tmax * length_of_L * dH) / ((length_of_L + K) * dS)                                                       # 356.17 K


        # dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_initial_den * dS)
        #
        # dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_initial_den * dS)
        #
        # dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_initial_den * dS)
        #
        # dGs[3] = (Tm_enzyme - T_initial_den) * dS
        #
        # dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (T_initial_den * dS)
        #
        # dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (T_initial_den * dS)
        #
        # dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (T_initial_den * dS)
        #
        # dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (T_initial_den * dS)
        #
        # dGs = np.clip(dGs, a_min=None, a_max=1e+12)
        #
        #
        # integration_initial_den = odeint(PCR_reaction, overall_concentration[e], time[0: t_initial_den * steps], args=(T_initial_den, dGs), mxstep=5000000)
        #
        # concentration[e, 0: t_initial_den * steps] = integration_initial_den
        #
        # overall_concentration[e] = integration_initial_den[-1]


        for i in range(number_cycles):

            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tden * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tden * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tden * dS)

            dGs[3] = (Tm_enzyme - Tden) * dS

            dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Tden * dS)

            dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Tden * dS)

            dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Tden * dS)

            dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Tden * dS)


            dGs = np.clip(dGs, a_min=None, a_max=1e+12)


            # integration_den = odeint(PCR_reaction_with_misbinding, overall_concentration[e], time[(t_initial_den * steps -1 + total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)], args=(Tden, dGs),  mxstep=5000000)
            #
            # concentration[e, t_initial_den * steps -1 + (total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)] = integration_den

            integration_den = odeint(PCR_reaction_with_misbinding, overall_concentration[e], time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

            concentration[e, (total * i * steps): ((total * i + tden) * steps)] = integration_den




            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tanneal * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tanneal * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tanneal * dS)

            dGs[3] = (Tm_enzyme - Tanneal) * dS

            dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Tanneal * dS)

            dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Tanneal * dS)

            dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Tanneal * dS)

            dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Tanneal * dS)


            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_anneal = odeint(PCR_reaction_with_misbinding, integration_den[-1], time[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs),  mxstep=5000000)
            #
            # concentration[e, t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)] = integration_anneal

            integration_anneal = odeint(PCR_reaction_with_misbinding, integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)

            concentration[e, ((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal


            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Text * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Text * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Text * dS)

            dGs[3] = (Tm_enzyme - Text) * dS

            dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Text * dS)

            dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Text * dS)

            dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Text * dS)

            dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Text * dS)


            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_ext = odeint(PCR_reaction_with_misbinding, integration_anneal[-1], time[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)], args=(Text, dGs),  mxstep=5000000)
            #
            # concentration[e, t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)] = integration_ext

            integration_ext = odeint(PCR_reaction_with_misbinding, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)
            concentration[e, ((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext



            overall_concentration[e] = integration_ext[-1]


        # dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_cooling_down * dS)
        #
        # dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_cooling_down  * dS)
        #
        # dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_cooling_down  * dS)
        #
        # dGs[3] = (Tm_enzyme - T_cooling_down ) * dS
        #
        # dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (T_cooling_down  * dS)
        #
        # dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (T_cooling_down  * dS)
        #
        # dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (T_cooling_down * dS)
        #
        # dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (T_cooling_down  * dS)
        #
        # dGs = np.clip(dGs, a_min=None, a_max=1e+12)
        #
        # integration_cool = odeint(PCR_reaction_with_misbinding, overall_concentration[e], time[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)], args=(T_cooling_down, dGs), mxstep=5000000)
        #
        # concentration[e, t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)] = integration_cool[-1]


        #values = integration_cool[-1]


    yield_PCR = ["S", "L"]

    concentration_PCR = np.zeros((len(overall_concentration), number_time_points, len(yield_PCR)))         # the matrix s length is all time points

    indexes_of_species = [[] for i in range(len(yield_PCR))]


    for x in range(len(yield_PCR)):


        for i in new_species:

            if yield_PCR[x] in i:


                indexes_of_species[x].append(new_species.index(i))


    print(indexes_of_species)



    for e in range(len(overall_concentration)):



        for x in range(number_time_points):


            for i in range(len(indexes_of_species)):


                summary_concentration = 0


                for n in range(len(indexes_of_species[i])):


                    summary_concentration = summary_concentration + concentration[e, x, indexes_of_species[i][n]]



                concentration_PCR[e, x, i] = summary_concentration



    yield_sum = np.zeros((len(overall_concentration), number_time_points))

    purity = np.zeros((len(overall_concentration), number_time_points))



    for e in range(len(overall_concentration)):


        for i in range(number_time_points):

            #print("S", concentration_PCR[e, i, 0])

            #print("L", concentration_PCR[e, i, 1])




            yield_sum[e, i] = concentration_PCR[e, i, 0] +  concentration_PCR[e, i, 1]

            purity[e, i] = concentration_PCR[e, i, 0] / yield_sum[e, i]



    for i in range(len(overall_concentration)):

        print("purity level of concentration set:", purity[i][-1])


    final_pur = []

    for i in range(len(overall_concentration)):

        final_pur.append(purity[i][-1])



    plt.figure(1)

    plt.suptitle("Purity level with different initial S1S2 values" , fontsize = 18, fontweight= 'bold')

    #plt.suptitle("Purity level with different annealing temperatures" , fontsize = 22, fontweight= 'bold')

    #plt.suptitle("Purity level with different lengths of primer mismatch\n after primer misbinding" , fontsize = 18, fontweight= 'bold')

    #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]



    for e in range(len(overall_concentration)):



        plt.plot(time, purity[e], label = label_s1s2[e])


    plt.legend(loc='upper left', prop={'size':8}, bbox_to_anchor=(1,1))

    plt.tick_params(labelsize = 12)

    plt.xlabel("Time (s) ",  FontSize= 16, fontweight= 'bold')
    plt.ylabel("Purity level: [S] / ( [S] + [L])",  FontSize= 16, fontweight= 'bold')


    # plt.figure(2)
    #
    # plt.suptitle("Purity level with different annealing temperatures" , fontsize = 22, fontweight= 'bold')
    #
    # plt.plot(T_multiple_annealing, final_pur,'ro')
    #
    # plt.tick_params(labelsize = 16)
    #
    # plt.xlabel("Annealing temperature (K) ",  FontSize= 18, fontweight= 'bold')
    # plt.ylabel("Purity level: [S] / ( [S] + [L]",  FontSize= 18, fontweight= 'bold')


    plt.show()




def purity_over_total_yield(all_concentration, time_vector):

    """ This function calculates the final purity of dsDNA throughout the PCR reaction
    primer misbinding"""


    yield_PCR = ["S", "L"]

    concentration_PCR = np.zeros((all_concentration.shape[0], len(yield_PCR)))         # the matrix s length is all time points

    indexes_of_species = [[] for i in range(len(yield_PCR))]


    for x in range(len(yield_PCR)):


        for i in new_species:

            if yield_PCR[x] in i:


                indexes_of_species[x].append(new_species.index(i))



    for x in range(all_concentration.shape[0]):


        for i in range(len(indexes_of_species)):

            summary_concentration = 0


            for n in range(len(indexes_of_species[i])):


                summary_concentration = summary_concentration + all_concentration[x, indexes_of_species[i][n]]



            concentration_PCR[x, i] = summary_concentration


    yield_sum = [0 for i in range(all_concentration.shape[0])]

    purity = [0 for i in range(all_concentration.shape[0])]




    for i in range(all_concentration.shape[0]):

        yield_sum[i] = concentration_PCR[i, 0] +  concentration_PCR[i, 1]

        purity[i] = concentration_PCR[i, 0] / yield_sum[i]


    print("Purity level at the end of PCR reaction with primer misbinding:", purity[-1])

    plt.figure(1)

    plt.suptitle("Purity level throughout the PCR\n with primer misbinding" , fontsize = 14, fontweight= 'bold')


    #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]


    #for i in range(len(yield_PCR)):

        #plt.subplot(2, 6, i+1)


        #plt.plot(time_vector, concentration_PCR[:, i], label = yield_PCR[i])      #label = misbinding_single_species_PCR[i]

        #plt.gca().set_title(misbinding_single_species_PCR[i])

        #plt.ylim([0, y_top_limit[i]])


    #plt.plot(time_vector, yield_sum, label = "Total yield: S + L")

    plt.plot(time_vector, purity)

    # plt.legend(loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))
    #
    # plt.xlabel("Time (s)")
    # plt.ylabel("Purity level: [S] / ( [S] + [L])")


    plt.tick_params(labelsize = 12)

    plt.xlabel("Time (s) ",  FontSize= 16, fontweight= 'bold')
    plt.ylabel("Purity level: [S] / ( [S] + [L])",  FontSize= 16, fontweight= 'bold')




    # plt.figure(2)
    #
    # plt.suptitle("Total concentrations of 4 single species over time with misbinding PCR" , fontsize = 14)
    #
    # highlighted_species = [0, 3, 5, 6]            # "S1", "S2", "P1", "P2", "Q1", "Q2", "E", "M2", "N2", "L1", "L2"
    #
    #
    # #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]
    #
    #
    # for i in range(len(highlighted_species)):
    #
    #     #plt.subplot(2, 4, i+1)
    #
    #
    #     plt.plot(time_vector, concentration_PCR[:, highlighted_species[i]], label = misbinding_single_species_PCR[highlighted_species[i]])
    #
    #     #plt.ylim([0, y_top_limit[i]])
    #
    #
    #     #plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})
    #
    # plt.xlabel("Time")
    # plt.ylabel("Total concentration")
    # plt.legend(loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))
    #
    #
    #
    # plt.figure(3)
    #
    # plt.suptitle("Total concentrations of misbinding single species over time" , fontsize = 14)
    #
    # highlighted_species = [0, 7, 8, 10]            # "S1", "S2", "P1", "P2", "Q1", "Q2", "E", "M2", "N2", "L1", "L2"
    #
    #
    # #y_top_limit = [6, 6, 8.3, 2.5, 0.22, 0.12, 10400, 0.11, 0.12, 0.1]
    #
    #
    # for i in range(len(highlighted_species)):
    #
    #     #plt.subplot(2, 4, i+1)
    #
    #
    #     plt.plot(time_vector, concentration_PCR[:, highlighted_species[i]], label = misbinding_single_species_PCR[highlighted_species[i]])
    #
    #     #plt.ylim([0, y_top_limit[i]])
    #
    #
    #     #plt.legend([species[plots[i]]], loc='upper left', prop={'size':10})
    #
    # plt.xlabel("Time")
    # plt.ylabel("Total concentration")
    # plt.legend(loc='upper left', prop={'size':11}, bbox_to_anchor=(1,1))

    plt.show()




def purity_diff_Tanneal(values):

    """ This function shows the impact of different annealing temperatures on final dsDNA purity """


    #Tanneal_scale = np.linspace(celsius_to_Kelvin(30), celsius_to_Kelvin(70), 41)

    Tanneal_scale = np.linspace(celsius_to_Kelvin(55), celsius_to_Kelvin(80), 26)

    dGs = [0 for x in range(8)]

    concentration = np.zeros((len(Tanneal_scale), number_time_points, len(new_species)))


    yield_PCR = ["S", "L"]

    concentration_PCR = np.zeros((len(Tanneal_scale), number_time_points, len(yield_PCR)))         # the matrix s length is all time points

    indexes_of_species = [[] for i in range(len(yield_PCR))]


    for x in range(len(yield_PCR)):


        for i in new_species:

            if yield_PCR[x] in i:


                indexes_of_species[x].append(new_species.index(i))


    #for e in range(len(overall_concentration)):

    for e in range(len(Tanneal_scale)):

        Tanneal = Tanneal_scale[e]

        values = initial_fixed

        # dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_initial_den * dS)
        #
        # dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_initial_den * dS)
        #
        # dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_initial_den * dS)
        #
        # dGs[3] = (Tm_enzyme - T_initial_den) * dS
        #
        # dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (T_initial_den * dS)
        #
        # dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (T_initial_den * dS)
        #
        # dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (T_initial_den * dS)
        #
        # dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (T_initial_den * dS)
        #
        # dGs = np.clip(dGs, a_min=None, a_max=1e+12)
        #
        #
        # integration_initial_den = odeint(PCR_reaction, overall_concentration[e], time[0: t_initial_den * steps], args=(T_initial_den, dGs), mxstep=5000000)
        #
        # concentration[e, 0: t_initial_den * steps] = integration_initial_den
        #
        # overall_concentration[e] = integration_initial_den[-1]


        for i in range(number_cycles):

            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tden * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tden * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tden * dS)

            dGs[3] = (Tm_enzyme - Tden) * dS

            dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Tden * dS)

            dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Tden * dS)

            dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Tden * dS)

            dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Tden * dS)


            dGs = np.clip(dGs, a_min=None, a_max=1e+12)


            # integration_den = odeint(PCR_reaction_with_misbinding, overall_concentration[e], time[(t_initial_den * steps -1 + total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)], args=(Tden, dGs),  mxstep=5000000)
            #
            # concentration[e, t_initial_den * steps -1 + (total * i * steps): t_initial_den * steps + ((total * i + tden) * steps)] = integration_den

            integration_den = odeint(PCR_reaction_with_misbinding, values, time[(total * i * steps): ((total * i + tden) * steps)], args=(Tden, dGs), mxstep=5000000)

            concentration[e, (total * i * steps): ((total * i + tden) * steps)] = integration_den




            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Tanneal * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Tanneal * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Tanneal * dS)

            dGs[3] = (Tm_enzyme - Tanneal) * dS

            dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Tanneal * dS)

            dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Tanneal * dS)

            dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Tanneal * dS)

            dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Tanneal * dS)


            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_anneal = odeint(PCR_reaction_with_misbinding, integration_den[-1], time[t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs),  mxstep=5000000)
            #
            # concentration[e, t_initial_den * steps + ((total * i + tden) * steps) - 1: t_initial_den * steps + ((total * i + tden + tanneal) * steps)] = integration_anneal

            integration_anneal = odeint(PCR_reaction_with_misbinding, integration_den[-1], time[((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)], args=(Tanneal, dGs), mxstep=5000000)

            concentration[e, ((total * i + tden) * steps) - 1: ((total * i + tden + tanneal) * steps)] = integration_anneal


            dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (Text * dS)

            dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (Text * dS)

            dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (Text * dS)

            dGs[3] = (Tm_enzyme - Text) * dS

            dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (Text * dS)

            dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (Text * dS)

            dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (Text * dS)

            dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (Text * dS)


            dGs = np.clip(dGs, a_min=None, a_max=1e+12)

            # integration_ext = odeint(PCR_reaction_with_misbinding, integration_anneal[-1], time[t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)], args=(Text, dGs),  mxstep=5000000)
            #
            # concentration[e, t_initial_den * steps + ((total * i + tden + tanneal) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps)] = integration_ext

            integration_ext = odeint(PCR_reaction_with_misbinding, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)], args=(Text, dGs), mxstep=5000000)
            concentration[e, ((total * i + tden + tanneal) * steps) - 1: (total * (i + 1) * steps)] = integration_ext



            values = integration_ext[-1]


        # dGs[0] = (Tmax * amplicon_length * dH) / ( amplicon_length + K) - (T_cooling_down * dS)
        #
        # dGs[1] = (Tmax * primer_length * dH) / ( primer_length + K) - (T_cooling_down  * dS)
        #
        # dGs[2] = (Tmax * extended_primer * dH) / ( extended_primer + K) - (T_cooling_down  * dS)
        #
        # dGs[3] = (Tm_enzyme - T_cooling_down ) * dS
        #
        # dGs[4] = (Tmax * length_misbinding_primer * dH) / ( length_misbinding_primer + K) - (T_cooling_down  * dS)
        #
        # dGs[5] = (Tmax * length_misbinding_extended_primer * dH) / ( length_misbinding_extended_primer + K) - (T_cooling_down  * dS)
        #
        # dGs[6] = (Tmax * length_misbinding_single_substrate * dH) / ( length_misbinding_single_substrate + K) - (T_cooling_down * dS)
        #
        # dGs[7] = (Tmax * length_of_L * dH) / ( length_of_L + K) - (T_cooling_down  * dS)
        #
        # dGs = np.clip(dGs, a_min=None, a_max=1e+12)
        #
        # integration_cool = odeint(PCR_reaction_with_misbinding, overall_concentration[e], time[t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)], args=(T_cooling_down, dGs), mxstep=5000000)
        #
        # concentration[e, t_initial_den * steps + (total * (i + 1) * steps) - 1: t_initial_den * steps + (total * (i + 1) * steps + t_cooling_down * steps)] = integration_cool[-1]


        #values = integration_cool[-1]

    for e in range(len(Tanneal_scale)):



        for x in range(number_time_points):


            for i in range(len(indexes_of_species)):


                summary_concentration = 0


                for n in range(len(indexes_of_species[i])):


                    summary_concentration = summary_concentration + concentration[e, x, indexes_of_species[i][n]]



                concentration_PCR[e, x, i] = summary_concentration



    yield_sum = np.zeros((len(Tanneal_scale), number_time_points))

    purity = np.zeros((len(Tanneal_scale), number_time_points))



    for e in range(len(Tanneal_scale)):


        for i in range(number_time_points):


            yield_sum[e, i] = concentration_PCR[e, i, 0] +  concentration_PCR[e, i, 1]

            purity[e, i] = concentration_PCR[e, i, 0] / yield_sum[e, i]

    final_purities = []

    for i in range(len(Tanneal_scale)):

        final_purities.append(purity[i][-1])


    final_yield_sum = []

    for i in range(len(Tanneal_scale)):

        final_yield_sum.append(yield_sum[i][-1])


    plt.figure(1)

    plt.suptitle("Purity level with different annealing temperatures" , fontsize = 22, fontweight= 'bold')

    plt.plot(Tanneal_scale, final_purities,'bo')

    plt.tick_params(labelsize = 16)

    plt.xlabel("Annealing temperature (K) ",  FontSize= 18, fontweight= 'bold')
    plt.ylabel("Purity level: [S] / ( [S] + [L])",  FontSize= 18, fontweight= 'bold')




    plt.figure(2)

    plt.suptitle("Product yield with different annealing temperatures" , fontsize = 22, fontweight= 'bold')

    plt.plot(Tanneal_scale, final_yield_sum,'ro')

    plt.tick_params(labelsize = 16)

    plt.xlabel("Annealing temperature (K) ",  FontSize= 18, fontweight= 'bold')
    plt.ylabel("Product yield (uM): [S] + [L]",  FontSize= 18, fontweight= 'bold')

    plt.show()



experiment_results = [121.01, 100.21, 70.7, 99.32, 30.19, 26.35]     # in ng/ul


new_model_result = [29.38990637941698, 3.984935187159642, 21.553138193134494, 17.60302588762814, 23.24906009323355, 25.24729240756307]      # in ng/ul

new_l1l2 = [1.4411534012276535, 0.04784522088128262, 0.8954339111076808, 0.5538091921659466, 1.7463379448408018, 3.1009165376022203]        # in ng/ul

new_purity = [0.9372456662910413, 0.9841191360943705, 0.9398215547270866, 0.9536681619545211, 0.9100037576361965, 0.8547099033283474]


# previous results

model_results = [13.399389300562468, 10.345881289354317, 30.30019959991723, 31.88182755512074, 3.9823027827329027, 7.696303301432735]       # in ng/ul

l1l2_results = [0.2580720786117887, 0.11973914167183851, 1.0865291921937494, 0.890811185488324, 0.14035884599348175, 0.511316613338765]     # in ng/ul

purity_results = [0.967225734185978, 0.9842826616815208, 0.9476291509779476, 0.9594607853861616, 0.9322468481479678, 0.8752663993704469]


def comparison_experiment_model():

    """ This function creates a plot with the results of laboratory experiments and the model's estimations """


    label_exp = ["Exp 0", "Exp 1", "Exp 2", "Exp 3", "Exp 4", "Exp 5"]


    plt.title("Comparison between dsDNA concentration obtained by\nlaboratory experiment and computational model ", FontSize= 22, FontWeight = "bold")        # position=(0.5, 1.05)

    plt.plot(range(len(experiment_results)), experiment_results,'bo',label ="Lab\nExperiment" )
    plt.plot(range(len(experiment_results)), experiment_results, 'b')
    plt.plot(range(len(model_results)), new_model_result, 'go', label = "Model")
    plt.plot(range(len(model_results)), new_model_result, 'g')

    plt.tick_params(labelsize = 16)

    plt.xticks(range(len(label_exp)), label_exp)

    plt.xlabel("Number of experiment", FontSize= 18, FontWeight = "bold")
    plt.ylabel("DsDNA concentration at\nthe end of PCR (ng/ul)", FontSize= 18, FontWeight = "bold")

    plt.legend(loc='upper left', prop={'size':16}, bbox_to_anchor=(1,1))

    plt.show()


if __name__ == '__main__':

    print("Second part of the code")



