from scipy.integrate import odeint

import matplotlib.pyplot as plt

import numpy as np



species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]


values = [0 for i in range(17)]

#values[0] = 3.0769e-15      # concentration of plasmid (S1S2) in mol/ul

values[0] = 1.18436669e-02
values[1] = 1.89253331e-02     # concentration of S1 in ng
values[2] = 1.89253331e-02        # concentration of S2 in ng
#values[3] = 3.0769e-10   # concentration of P1 in mol/ul
#values[4] = 3.0769e-10   # concentration of P2 in mol/ul

values[3] = 0.2
values[4] = 0.2


values[5] = 0       # concentration of S1P2
values[6] = 0       # concentration of S2P1
values[7] = 200  # concentration of E in Unit
values[8] = 0       # concentration of S1P2E
values[9] = 0       # concentration of S2P1E
values[10] = 200 # concentration of dNTP in mol/ul
values[11] = 0      # concentration of Q1
values[12] = 0      # concentration of Q2
values[13] = 0      # concentration of S1Q2E
values[14] = 0      # concentration of S2Q1E
values[15] = 0      # concentration of S2Q1
values[16] = 0      # concentration of S1Q2

Tden = 310          # Kelvin which is equal to 90 degree Celsius

#Tden = 369.15  # 96 degree

#Tanneal = 345.15        # Kelvin which is equal to 72 degree Celsius

Tanneal = 303.15        # 30 degree
#Tanneal = 369.15

Text = 303.15           # Kelvin which is equal to 80 degree Celsius

#Text = 343.15  # 70 degree

#Text = 369.15

tden = 10           # seconds
tanneal = 10      # seconds
text = 10            # seconds


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



R = 8.314e-3        # Gas contant in  J / K mol

#kB = 1.38064852e-23  # Boltzmann constant"


#

#print("The number of added nucleotides to the primer during primer_ext_1:", n)



Tm_S1S2 = 358.15

Tm_primer = 308.15


dS = -1256.0395 # j/m

#dS = - 52.7184


def denaturation(values, t, T, dGs):


    S1S2 = values[0]
    S1 = values[1]
    S2 = values[2]

    #print("S1S2, S1, S2", S1S2, S1, S2 )


    kf1 = 1



    exponent_1 = dGs[0]/(R*T)

    exponent_1 = np.clip(exponent_1, a_min= None,  a_max= 29.9336 )

    kr1 = kf1 * np.exp(exponent_1)

    #print("kr1 origin", kr1)

    kr1 = np.clip(kr1, a_min= -1e+14, a_max= 1e+14)

    #print("kr1 second", kr1)

    #kr1 = 10




    rate_den = - kr1 * S1S2 + kf1 * S1 * S2

    #print("rateden first", rate_den)

    rate_den = np.clip(rate_den,  a_min= -1e+14, a_max= 1e+14 )

    #print("rateden", rate_den)

    y = np.zeros(17)

    y[0] = rate_den
    #y[0] = np.clip(y[0],  a_min= 0, a_max= 1e+14 )
    y[1] = -rate_den
    #y[1] = np.clip(y[1],  a_min= 0, a_max= 1e+14 )
    y[2] = -rate_den
    #y[2] = np.clip(y[2],  a_min= 0, a_max= 1e+14 )


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

    #print("exponent_2a origin", exponent_2a)

    exponent_2a = np.clip(exponent_2a, a_min= None, a_max= 29.9336)

    #print("exponent_2a second", exponent_2a)

    kr2a = kf2 * np.exp(exponent_2a)

    kr2a = np.clip(kr2a,  a_min= -1e+14, a_max= 1e+10 )



    print("kr2a second", kr2a)



    rate_S1P2_bind = kf2 * S1 * P2 - kr2a * S1P2


    rate_S1P2_bind = np.clip(rate_S1P2_bind,  a_min= -1e+14, a_max= 1e+10 )

    if np.abs(rate_S1P2_bind<1e-14):
        rate_S1P2_bind = 0

    print("rate S1P1 bind", rate_S1P2_bind)



    exponent_2b = dGs[1]/(R*T)

    #print("exponent_2b origin", exponent_2b)

    exponent_2b = np.clip(exponent_2b, a_min= None, a_max= 29.9336)

    #print("exponent_2b second", exponent_2b)

    kr2b = kf2 * np.exp(exponent_2b)


    kr2b = np.clip(kr2b,  a_min= -1e+14, a_max= 1e+10 )

    #print("kr2b", kr2b)

    rate_S2P1_bind = kf2 * S2 * P1 - kr2b * S2P1


    rate_S2P1_bind = np.clip(rate_S2P1_bind,  a_min= -1e+14, a_max= 1e+10 )

    if np.abs(rate_S2P1_bind<1e-14):
        rate_S2P1_bind = 0

    print("rate_s2P1", rate_S2P1_bind)


    y = np.zeros(17)


    y[1] = - rate_S1P2_bind
    y[2] = - rate_S2P1_bind
    y[3] = - rate_S2P1_bind
    y[4] = - rate_S1P2_bind
    y[5] = rate_S1P2_bind
    y[6] = rate_S2P1_bind


    return y


def first(values, t, T, dGs):

    summary = denaturation(values, t, T, dGs) + primer_binding_1(values, t, T, dGs)

    return summary



# time_first = np.linspace(0, tden, tden +1)
#
# dGs = [0,0,0,0]
#
# dGs[0] = (Tm_S1S2 - Tden) * dS
#
# print(dGs)
#
#
# integration_first = odeint(denaturation, values, time_first , args=(Tden, dGs))
#
#
#
#
#
#
#
#     #print(denaturation(values, 10, 340, [-300, 0,0,0]))
#
#
#     #integration_first = odeint(denaturation, values, time_first , args=(Tden, dGs))
#
# print(integration_first)
#
# plt.figure(1)
#
# plt.plot(time_first, integration_first[:, 0])
# plt.plot(time_first, integration_first[:, 1])
# plt.plot(time_first, integration_first[:, 2])


# integration_sec = odeint(primer_binding_1, values, time_first , args=(Tden, dGs))
#     #
# print(integration_sec)
#     #
# plt.figure(2)
#     #
# plt.plot(time_first, integration_sec[:, 1])
# plt.plot(time_first, integration_sec[:, 3])
# plt.plot(time_first, integration_sec[:, 5])
#
# plt.show()
#
# integration = odeint(first, values, time_first , args=(Tden, dGs))
#
# plt.figure(1)
#
# plt.plot(time_first, integration[:, 0])
# plt.plot(time_first, integration[:, 1])
# plt.plot(time_first, integration[:, 2])
#
# plt.figure(2)
#
# plt.plot(time_first, integration[:, 1])
# plt.plot(time_first, integration[:, 3])
# plt.plot(time_first, integration[:, 5])
#
# plt.show()


dGs = [0,0,0,0]

concentration = np.empty((number_time_points, 17))

for i in range(number_cycles):




    dGs[0] = (Tm_S1S2 - Tden) * dS

    dGs[1] = (Tm_primer - Tden) * dS


    dGs = np.clip(dGs, a_min= None, a_max= 1e+14 )

    print("dGs_den", dGs)



    integration_den = odeint(primer_binding_1, values, time[(total * i * steps) : ((total * i + tden) * steps)] , args=(Tden, dGs))

    concentration[(total * i * steps) : ((total * i + tden) * steps)] = integration_den





    dGs[0] = (Tm_S1S2 - Tanneal) * dS

    dGs[1] = (Tm_primer - Tanneal) * dS

    dGs = np.clip(dGs, a_min= None, a_max= 1e+14 )

    print("dGs_anneals", dGs)

    integration_anneal = odeint(primer_binding_1, integration_den[-1], time[((total * i + tden) * steps)-1 : ((total * i + tden + tanneal) * steps)] , args=(Tanneal, dGs ))

    concentration[((total * i + tden) * steps) - 1 : ((total * i + tden + tanneal) * steps)] = integration_anneal


    dGs[0] = (Tm_S1S2 - Text) * dS

    dGs[1] = (Tm_primer - Text) * dS

    dGs = np.clip(dGs, a_min= None, a_max= 1e+14 )

    print("dGs_text", dGs)

    integration_ext = odeint(primer_binding_1, integration_anneal[-1], time[((total * i + tden + tanneal) * steps) -1 : (total * (i + 1) * steps)], args=(Text, dGs))

    concentration[((total * i + tden + tanneal) * steps) -1 : (total * (i + 1) * steps)] = integration_ext


    values = integration_ext[-1]

    print(values)


    #print(concentration)


print("whyyy",(Tm_S1S2 - Tden, dS, (Tm_S1S2 - Tden) * dS))

plt.plot(time, concentration[:, 0])


#plt.legend(["S1S2"], loc='lower left')

plt.plot(time, concentration[:, 1])

plt.plot(time, concentration[:, 3])

plt.plot(time, concentration[:, 5])



plt.legend(["S1S2", "S1", "primer", "SP"], loc='best')




plt.show()



#
#     plt.figure(1)
#
#     for i in range(8):
#
#         plt.subplot(2, 4, i+1)
#
#
#         plt.plot(time, concentration[:, i])
#
#         plt.legend([species[i]], loc='best', prop={'size':10})
#
#         plt.xlabel("Time")
#         plt.ylabel("Concentration")
#
#
#     plt.figure(2)
#
#     plt.suptitle("Change of the species' concentrations over time" , fontsize = 14)
#
#     for i in range(9):
#
#
#         plt.subplot(2, 5, i +1)
#
#
#         plt.plot(time, concentration[:, i + 8])
#
#         plt.legend([species[i + 8]], loc='best', prop={'size':10})
#
#         plt.xlabel("Time")
#         plt.ylabel("Concentration")

plt.show()

