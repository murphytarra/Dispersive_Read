
import ECHO 
from ECHO import ECHO_3x3 as e3
from ECHO import ECHO_2x2 as e2
Q = e3.ECHO_3x3()
Q.print_parameters()
Q.run()
Q.coher_sweep(B = 10, gmin = 0.25, gmax = 2, emin = -4, emax = 4)
Q.plot_eigen()
Q.freq_sweep()
Q.coher_sweep()

P = e2.ECHO_2x2()
P.run()
