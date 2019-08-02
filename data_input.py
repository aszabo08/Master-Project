


import numpy as np



# ----- Functions for converting units --------------------------------------------------------------------------------------------------------------------


# from previous code --- modified: from nanomolar to micromolar


def mass_to_micromolar(ng, bp, V_ul) :
	"""
	The nanomolar concentration resulting when ng nanograms of dsDNA (linear or circular)
	of length bp base pairs is put in microlitre volume V_ul
	"""

	# 1 mol of single base pairs weights 660 grams
	mol = (ng / 1e9) / (660 * bp)	# grams / grams per mole of dsDNA molecules = number of moles of molecules
	M = mol / (V_ul / 1e6)
	return M * 1e6


#### modified!!! from nano to micro

def U_to_micromolar(U) :
	"""
	Guestimate:
	We assume 1.25U per 50uL to correspond to the nM amount of enzyme
	which gives an acceptable amplification curve over 35 cycles,
	which is about 200nM
	"""
	return U * (0.2 / 1.25)




def uM_primer(ul):

    "Assuming that the stock concentration is 10.0 uM and the final volume is 50 ul"


    return (ul * 10) / 50


def uM_dNTP(ul):

    "Assuming that the stock concentration is 10.0 mM (10000 uM) and the final volume is 50 ul"


    return (ul * 10000) / 50



def celsius_to_Kelvin(x):

    return x + 273.15




# ----- Constant values ------------------------------------------------------------------------------------------------------------------------------------



# PART I : PCR with accurate primer binding site



species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]


values = [0 for i in range(33)]



values[0] = float(input("Enter the concentration of plasmid (ng): "))

values[0] = mass_to_micromolar(values[0], 1000, 50)


values[3] = float(input("Enter the concentration of each primer (uL): "))

values[3] = uM_primer(values[3])

values[4] = values[3]


values[7] = float(input("Enter the concentration of polymerase (U): "))

values[7] = U_to_micromolar(values[7])


values[10] = float(input("Enter the concentration of each dNTP (uL): "))

values[10] = uM_dNTP(values[10])



Tden_celsius, tden_string = input("Enter the temperature of denaturation (Celsius) and the length of it (second) ").split()

# Converting the temperature from Celsius to Kelvin

Tden = float(Tden_celsius)

Tden = celsius_to_Kelvin(Tden)

tden = int(tden_string)


Tanneal_celsius, tanneal_string = input("Enter the temperature of annealing (Celsius) and the length of it (second) ").split()

# Converting the temperature from Celsius to Kelvin

Tanneal = float(Tanneal_celsius)

Tanneal = celsius_to_Kelvin(Tanneal)

tanneal = int(tanneal_string)



Text_celsius, text_string = input("Enter the temperature of primer extension (Celsius) and the length of it (second) ").split()

# Converting the temperature from Celsius to Kelvin

Text = float(Text_celsius)

Text = celsius_to_Kelvin(Text)

text = int(text_string)

number_cycles = int(input("Enter the number of cycles "))



T_cooling_down = 298.15             # 25 degree celsius

t_cooling_down = 10





initial_dNTP = values[10]

initial_fixed = values

functions_name = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer_ext_2"]



total = tden + tanneal + text



steps = 1


amplicon_length = 1000

primer_length = 15


#primer_length = 27



n = 10

extended_primer = primer_length + n

extended_length = amplicon_length - extended_primer

time = np.linspace(0, number_cycles * (tden + tanneal + text) + t_cooling_down, number_cycles * (tden + tanneal + text) * steps + t_cooling_down * steps)  # for every second "steps" points are distinguished

number_time_points = total * steps * number_cycles + t_cooling_down * steps


dS = - 2


max_exponent = 15    # 13

min_clip = -1e+15       #18

max_clip = 1e+15



forward_rate = 1


R = 8.314e-3        # Gas contant in  J / K mol


#Tm_primer = 320.15      # 47 deggree

#Tm_extended_primer = 337.15     #64 degree


Tm_primer = 343.15

Tm_extended_primer = 349.15


Tmax = 373.15           # 100 degree


K = (primer_length * extended_primer * (Tm_extended_primer - Tm_primer)) / (extended_primer * Tm_primer - primer_length * Tm_extended_primer)


dH = (Tm_primer * dS * (primer_length + K)) / (Tmax * primer_length)


dH_check = (Tm_extended_primer * dS * ( extended_primer + K)) / ( Tmax * extended_primer)

Tm_S1S2 = (Tmax * amplicon_length * dH) / ((amplicon_length + K) * dS)
#
# #Tm_enzyme = dH / dS


Tm_enzyme = 353.15              # 80 degree  not calculated!




# PART II : PCR with primer misbinding



new_species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1", "S1M2", "S1M2E", "S1N2E", "S1L2", "L2", "S1N2", "L2P1", "L2P1E", "L2Q1E", "L1L2", "L2Q1", "L1", "L1P2", "L1P2E", "L1Q2E", "L1Q2"]

# amplicon length - primer length - n = 1000 - 15 -10 = 975         -------->  the longest extended length is 975 -1, when the product is only one nucleotide short
#                                                                              the shortest extended length is 0
# average extended length = (974 + 0)/2 = 487

misbinding_extended_length = 487

length_of_L = 487 + 15 + 10         # 512


# mis = 0.25
#
#
# # Assuming that one-fourth of the full length is mismatched: we just deduct that length
#
# length_misbinding_primer = primer_length - round(mis * primer_length)
#
# length_misbinding_extended_primer = extended_primer - round(mis * extended_primer)
#
# length_misbinding_single_substrate = length_of_L - round(mis * length_of_L)


mismatch = 5


#Assuming that we have 5 mismatches in the primer: the extensions will be based on correct base pairing!

length_misbinding_primer = primer_length - mismatch

length_misbinding_extended_primer = extended_primer - mismatch

length_misbinding_single_substrate = length_of_L - mismatch


Tm_misbinding_primer = (Tmax * length_misbinding_primer * dH) / ((length_misbinding_primer + K) * dS)                                       # 306.1179269849045 K

Tm_misbinding_extended_primer = (Tmax * length_misbinding_extended_primer * dH) / ((length_misbinding_extended_primer + K) * dS)            # 328.8778398293736 K

Tm_misbinding_single_substrate = (Tmax * length_misbinding_single_substrate * dH) / ((length_misbinding_single_substrate + K) * dS)          # 364.27564234055836 K

Tm_misbinding_double_substrate = (Tmax * length_of_L * dH) / ((length_of_L + K) * dS)                                                        # 364.78659575718 K



















