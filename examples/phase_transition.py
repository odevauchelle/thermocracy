from pylab import *

import sys
sys.path.append('./../')

from thermocracy import population

pop = population( connectivity = ( rand(*[50]*2) > .5 )*1, epsilon = 1, beta = 2., state = None )

X_mean = [ mean( pop.get_state_vector() ) ]

for beta in [.5,1,2] :

    pop.beta = beta

    for _ in range(30) :

        pop.new_deal()

        for _ in range(100)   :
            pop.evolve(10)
            X_mean += [ mean( pop.get_state_vector() ) ]

    hist( X_mean, bins = linspace(-1,1,50), alpha = 0.5, label = str( beta ), normed = True )

show()
