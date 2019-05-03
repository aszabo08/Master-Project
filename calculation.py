




def creating_amplicon(amplicon_length):



    bases = ["A", "G", "C", "T"]

    amplicon_list = []

    for i in range(amplicon_length):

        amplicon_list.append(random.choice(bases))


    amplicon = "".join(amplicon_list)


    return amplicon




def creating_primer(primer_lenght):

    bases = ["A", "G", "C", "T"]

    primer_list = []

    for i in range(primer_lenght):

        primer_list.append(random.choice(bases))


    primer = "".join(primer_list)


    return primer






amplicon_length = 1000

    primer_length = 15

    amplicon = creating_amplicon(amplicon_length)   # the resulted amplicon: TACAGGGATAGTTCAACTTTAGGGCTCAAGCGGAAGTTGGCATCACCCAAGCGACTGGGGCAAATGATGAGGGGGAGTCCTTCGGGTTGATTACCTAACGTCTGTCTGTATCAGGCCCGCAAGCTATGTCTCCCTTCCGTGATGCATGGAGAACCTGCGCTCAAGGGGAAATTAGCTCGTACTCTTCGCGCAGGGGGCATGCTGTGCGGACTCAATTAGTTGTTTACTGGCTTGAGGAATTTTCGTCCGGTATATAATCACATGCAGTAAAACCCTGATAGCGGTTACTTCTTGAGCAAACTTTTACGTGTTCTTCGCAGGTACACAACTTCGCTACCTTGCATAGGCATGTGTATGCTGAAAGGACCTATGCACTAACATAACTTAGTAGTAGTGAGACAACTCGAATTCAAGCTATTCCTGCTGCAAAGAGATCACCTATCGTCGGTCTCCGAGGGCGTAAAGCCATCGAGATTACCAGACTGGTGGGCGATTCCATCGACGACGTCAGCCTTCAGACATTCTAATAGGACCTCTGGGGCTGACAATGAGAGGTCCTGTTCTGGATTTGTAAGAGCCTCATTGTGTCAGAACCACAATTGATATGATCGGTTTTAACTACAATCGGATCCACCAAAACTCCATGCTAGAGCCAAGGATAGCTCGGATGAAGTGTGTAAATCAGATACAACCCTTTCCTATAATCCTACGATATATACCGTGACATCGGGTGGCTCTCTCCCACCCCCGGCAGTAGACCAAGCAGTCCATCCCACTGAGCCATTGTGACATAGCTTGTAAGTATCATTCACTATAACGCAACGCCGGGTAGCCTCTACGGTCGTCCTGACTAGTACATAATTGTGGACCTCCATGAGGAGTACAGTGTCAACTTACTAGTCCCTGACTGTTCCGAACGTGTGCCTAAATTAAGACTGGAGCGAATATCCCCTGTCTCACAGTGAGACCACAACTAAAAGCGAGTCGTCCTACGTATGAG

                                                    # complement of amplicon: ATGTCCCTATCAAGTTGAAATCCCGAGTTCGCCTTCAACCGTAGTGGGTTCGCTGACCCCGTTTACTACTCCCCCTCAGGAAGCCCAACTAATGGATTGCAGACAGACATAGTCCGGGCGTTCGATACAGAGGGAAGGCACTACGTACCTCTTGGACGCGAGTTCCCCTTTAATCGAGCATGAGAAGCGCGTCCCCCGTACGACACGCCTGAGTTAATCAACAAATGACCGAACTCCTTAAAAGCAGGCCATATATTAGTGTACGTCATTTTGGGACTATCGCCAATGAAGAACTCGTTTGAAAATGCACAAGAAGCGTCCATGTGTTGAAGCGATGGAACGTATCCGTACACATACGACTTTCCTGGATACGTGATTGTATTGAATCATCATCACTCTGTTGAGCTTAAGTTCGATAAGGACGACGTTTCTCTAGTGGATAGCAGCCAGAGGCTCCCGCATTTCGGTAGCTCTAATGGTCTGACCACCCGCTAAGGTAGCTGCTGCAGTCGGAAGTCTGTAAGATTATCCTGGAGACCCCGACTGTTACTCTCCAGGACAAGACCTAAACATTCTCGGAGTAACACAGTCTTGGTGTTAACTATACTAGCCAAAATTGATGTTAGCCTAGGTGGTTTTGAGGTACGATCTCGGTTCCTATCGAGCCTACTTCACACATTTAGTCTATGTTGGGAAAGGATATTAGGATGCTATATATGGCACTGTAGCCCACCGAGAGAGGGTGGGGGCCGTCATCTGGTTCGTCAGGTAGGGTGACTCGGTAACACTGTATCGAACATTCATAGTAAGTGATATTGCGTTGCGGCCCATCGGAGATGCCAGCAGGACTGATCATGTATTAACACCTGGAGGTACTCCTCATGTCACAGTTGAATGATCAGGGACTGACAAGGCTTGCACACGGATTTAATTCTGACCTCGCTTATAGGGGACAGAGTGTCACTCTGGTGTTGATTTTCGCTCAGCAGGATGCATACTC

    primer = creating_primer(primer_length)         # the resulted primer: GAATGGTCGCTCGCG

                                                    # complement of the primer: CTTACCAGCGAGCGC



    values = [0 for i in range(12)]                 # Initializing the concentrations of the reactants at 0



    values[0] = float(input("Enter the concentration of plasmid (ng): "))


    values[3] = float(input("Enter the concentration of each primer (microL): "))


    values[4] = values[3]


    values[7] = float(input("Enter the concentration of polymerase (U): "))


    values[10] = float(input("Enter the concentration of each dNTP (microL): "))



    Tden_celsius, tden_string = input("Enter the temperature of denaturation (Celsius) and the length of it (second) ").split()

    # Converting the temperature from Celsius to Kelvin

    Tden = float(Tden_celsius) + 273.15

    tden = float(tden_string)


    Tanneal_celsius, tanneal_string = input("Enter the temperature of annealing (Celsius) and the length of it (second) ").split()

    # Converting the temperature from Celsius to Kelvin

    Tanneal = float(Tanneal_celsius) + 273.15

    tanneal = float(tanneal_string)



    Text_celsius, text_string = input("Enter the temperature of primer extension (Celsius) and the length of it (second) ").split()

    # Converting the temperature from Celsius to Kelvin

    Text = float(Text_celsius) + 273.15

    text = float(text_string)


    number_cycles = int(input("Enter the number of cycles "))
