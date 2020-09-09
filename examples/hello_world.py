from pylab import *

import sys
sys.path.append('./../')

from thermocracy import population

pop = population( connectivity = ( rand(*[50]*2) > .5 )*1, epsilon = 1, beta = 1, state = None )

X_mean = [ mean( pop.get_state_vector() ) ]

for _ in range(1000)   :
    pop.evolve(10)
    X_mean += [ mean( pop.get_state_vector() ) ]

hist( X_mean, bins = linspace(-1,1,20) )

show()
