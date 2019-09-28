from data_input import *
from PCR_project import *
from misbinding_primer import *

if __name__ == '__main__':

    if (t_initial_den > 0) and (t_cooling_down > 0):
        extended_PCR_misbinding_integration(values)
    else:
        PCR_misbinding_integration(values)
