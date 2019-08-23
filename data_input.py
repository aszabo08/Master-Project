

import numpy as np




# this piece of source code is originated from pcrparam.py, created by Dr Benjamin Shirt - Ediss
# I modified it to micromolar, instead of the original nanomolar concentration
def mass_to_micromolar(ng, bp, V_ul) :
	"""
	The micromolar concentration resulting when ng nanograms of dsDNA (linear or circular)
	of length bp base pairs is put in microlitre volume V_ul
	"""

	# 1 mol of single base pairs weights 660 grams
	mol = (ng / 1e9) / (660 * bp)	# grams / grams per mole of dsDNA molecules = number of moles of molecules
	M = mol / (V_ul / 1e6)
	return M * 1e6


# this piece of source code is originated from pcrparam.py, created by Dr Benjamin Shirt - Ediss
# I modified it to micromolar, instead of the original nanomolar concentration
def U_to_micromolar(U) :
	"""
	Guestimate:
	We assume 1.25U per 50uL to correspond to the nM amount of enzyme
	which gives an acceptable amplification curve over 35 cycles,
	which is about 200nM
	"""
	return U * (0.2 / 1.25)


def uM_to_ng_per_ul(uM, bp) :

	""" This function converts the micromolar concentration into ng/ul"""

	g_per_ul = uM * 1e-6 * 1e-6 * 660 * bp		# uM to M to mol/ul to g/ul

	return g_per_ul * 1e9




# def uM_primer(ul):
#
# 	""" Converting ul to uM, assuming that the stock concentration is 10.0 uM and the final volume is 50 ul """
#
# 	return (ul * 10) / 50
#
# def uM_dNTP(ul):
#
# 	""" Converting ul to uM, assuming that the stock concentration is 10.0 mM (10000 uM) and the final volume is 50 ul"""
#
#     return (ul * 10000) / 50



def celsius_to_Kelvin(x):

	""" Converting degree Celisus into Kelvin"""

	return x + 273.15




species = ["S1S2", "S1", "S2", "P1", "P2", "S1P2", "S2P1", "E", "S1P2E", "S2P1E", "dNTP", "Q1", "Q2", "S1Q2E", "S2Q1E", "S1Q2", "S2Q1"]


values = [0 for i in range(33)]



values[0] = float(input("Enter the concentration of plasmid (ng): "))

values[0] = mass_to_micromolar(values[0], 1000, 50)


#values[3] = float(input("Enter the concentration of each primer (uL): "))
#values[3] = uM_primer(values[3])

values[3] = float(input("Enter the concentration of each primer (uM): "))
values[4] = values[3]


values[7] = float(input("Enter the concentration of polymerase (U): "))
values[7] = U_to_micromolar(values[7])


#values[10] = float(input("Enter the concentration of each dNTP (uL): "))
#values[10] = uM_dNTP(values[10])

values[10] = float(input("Enter the concentration of each dNTP (uM): "))


T_initial_den_celsius, t_initial_den_string = input("Enter the temperature of initial denaturation (Celsius) and the length of it (second) ").split()

# Converting the temperature from Celsius to Kelvin

T_initial_den = float(T_initial_den_celsius)

T_initial_den = celsius_to_Kelvin(T_initial_den)

t_initial_den = int(t_initial_den_string)




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


T_cooling_down_celsius, t_cooling_down_string = input("Enter the temperature of final elongation (Celsius) and the length of it (second) ").split()

# Converting the temperature from Celsius to Kelvin

T_cooling_down = float(T_cooling_down_celsius)

T_cooling_down = celsius_to_Kelvin(T_cooling_down)

t_cooling_down = int(t_cooling_down_string)




number_cycles = int(input("Enter the number of cycles "))

enzyme_type = input("Enter the type of enyzme: 'taq' or 'Q5' ")


amplicon_length = 2207                     # bp
primer_length = 25                         # nt
n = 10                                     # nt

Tm_primer = 339.15                       # melting temperature of the primer in K (66 °C)
Tm_extended_primer = 346.15              # melting temperature of the extended primer in K (73 °C)

functions_name = ["denaturation", "primer_binding_1", "polymerase_binding_1", "primer_ext_1", "polymerase_binding_2", "primer_binding_2", "primer_ext_2"]

steps = 1


initial_fixed = values
total = tden + tanneal + text
extended_primer = primer_length + n
extended_length = amplicon_length - extended_primer
time = np.linspace(0, t_initial_den + number_cycles * (tden + tanneal + text) + t_cooling_down, number_cycles * (tden + tanneal + text) * steps + (t_cooling_down + t_initial_den) * steps)  # for every second "steps" points are distinguished
number_time_points = total * steps * number_cycles + (t_initial_den + t_cooling_down) * steps

dS = - 2                                    # entropy change, unit: kJ K^-1

max_exponent = 7
min_clip = -1e+12
max_clip = 1e+12

forward_rate = 1                                # uM/s
R = 8.314e-3                                    # Gas constant in kJ K^-1 mol^-1

Tmax = 373.15                                   # maximum melting temperature in K (100 °C)

K = (primer_length * extended_primer * (Tm_extended_primer - Tm_primer)) / (extended_primer * Tm_primer - primer_length * Tm_extended_primer)           # Michaelis constant

dH = (Tm_primer * dS * (primer_length + K)) / (Tmax * primer_length)                                 # enthalpy change, unit: kJ

Tm_S1S2 = (Tmax * amplicon_length * dH) / ((amplicon_length + K) * dS)

Tm_enzyme = 356.15                      # 85 °C

misbinding_extended_length = 487


# for the experiment
# misbinding_extended_length = 1086


#length_of_L = 487 + 15 + 10         # 512

length_of_L = misbinding_extended_length + primer_length + n

#mismatch = int(primer_length * 0.75)

mismatch = int(0.8 * primer_length)

#mismatch = int(0.8 * primer_length)

length_misbinding_primer = primer_length - mismatch

length_misbinding_extended_primer = extended_primer - mismatch

length_misbinding_single_substrate = length_of_L - mismatch



Tm_misbinding_primer = (Tmax * length_misbinding_primer * dH) / ((length_misbinding_primer + K) * dS)                                       # 342.71 K

Tm_misbinding_extended_primer = (Tmax * length_misbinding_extended_primer * dH) / ((length_misbinding_extended_primer + K) * dS)            # 346.87 K

Tm_misbinding_single_substrate = (Tmax * length_misbinding_single_substrate * dH) / ((length_misbinding_single_substrate + K) * dS)         # 356.168 K

Tm_misbinding_double_substrate = (Tmax * length_of_L * dH) / ((length_of_L + K) * dS)

