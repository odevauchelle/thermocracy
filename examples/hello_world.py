from pylab import *

import sys
sys.path.append('./../')

from thermocracy import population, plot_connectivity

pop = population( connectivity = ( rand(*[15]*2) > .5 )*1, epsilon = 1, beta = 2, state = None )

X_mean = [ mean( pop.get_state_vector() ) ]

for _ in range(100)   :
    pop.evolve(10)
    X_mean += [ mean( pop.get_state_vector() ) ]

figure()
plot_connectivity(pop.connectivity)

figure()
hist( X_mean, bins = linspace(-1,1,10) )

show()
