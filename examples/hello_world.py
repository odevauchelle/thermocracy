from pylab import *

import sys
sys.path.append('./../')

import thermocracy as thm

H = thm.Hamiltonian( terms = [ thm.neighbors_influence, thm.polls_influence ], coeffs = [ 1, 1 ] )

pop = thm.population( connectivity = ( rand(*[15]*2) > .5 )*1, H = H, beta = 2, state = None )

X_mean = [ mean( pop.get_state_vector() ) ]

for _ in range(100)   :
    pop.evolve(10)
    X_mean += [ mean( pop.get_state_vector() ) ]

figure()
thm.plot_connectivity(pop.connectivity)

figure()
hist( X_mean, bins = linspace(-1,1,10) )

show()
