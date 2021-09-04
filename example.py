
"""
Example code for ECHO 

3 Level System with all parameters set to initial values

For further information, see pdf 'Dispersive Readout of Multilevel Systems'

by Tara Murphy

Completed during 6 Week Summer Internship at Quantum Motion Technologies, UCL
"""

#%reset
from Dispersive_Readout_Package import ECHO_3x3 as dr
from Dispersive_Readout_Package import ECHO_2x2 as ds

Q = dr.ECHO_3x3()
Q.print_parameters()
Q.run()
Q.coher_sweep(B = 10, gmin = 0.25, gmax = 2, emin = -4, emax = 4)
Q.plot_eigen()
