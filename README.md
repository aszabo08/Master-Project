This repository contains the improving versions of my master's project for MSc Bioinformatics at Newcastle University from where I graduated in 2019.

The aim of this dissertation was to develop a computational model of polymerase chain reaction (PCR) including the possibility of primer misbinding, providing a more sensitive PCR model to support the teaching process of complex experimental design.

The created model is based on a dynamical system, displaying 33 species and 22 interactions between them. Primer misbinding can only occur on one of the DNA strands in order to keep the number of participating species relatively low. For the same reason only one extended primer stage was distinguished. The model includes certain thermodynamic parameters such as change of entropy, enthalpy and Gibbs free energy. The dNTP incorporation rate of the enzyme is set to be temperature-dependent, and the enzyme denatures at high temperatures. The set value of entropy difference, rate of enzyme denaturation and the effective length of primer mismatch was estimated after assessing different values’ impact on the behaviour of the model.

In order to run the PCR model, the user provides 16 experimental factors: the concentrations of the dsDNA, primers, enzyme and dNTP, the temperature and length of initial denaturation, denaturation, annealing, extension and final extension. Furthermore, the number of cycles and the type of enzyme is also chosen by the user. The enzyme type can be either Taq or Q5, which show different enzyme activity at different temperatures. The outcome is five plots displaying the concentration change over time of all the participant species and the purity level of the dsDNA. 


Instruction
-------------------------------------------------------------------------------------------
To run the program: download the Master-Project and run the "master_PCR_complete.py" file.
-------------------------------------------------------------------------------------------

An example input for the species' contentration and temperature setting in the polymerase chain reaction model:

plasmid (ng): 5

primer (uM): 0.5  

polymerse (U): 0.2

dNTP (uM): 200 

initial denaturation (°C and s): 96 30

denaturation (°C and s): 96 10

annealing (°C and s): 61 10

extension (°C and s): 72 20

final elongation (°C and s): 72 300

cumber of cycles: 32

enzyme type: taq
