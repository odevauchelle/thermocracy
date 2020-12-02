from pylab import *

import sys
sys.path.append('./../')

import thermocracy as thm

N = 50

H = thm.Hamiltonian( terms = [ thm.neighbors_influence, thm.polls_influence], coeffs = [ 1, 1 ] )

pop = thm.population( connectivity = ( rand(*[N]*2) > .5 )*1, H = H, beta = None, state = None )



for beta in array([.5])/N :

    pop.beta = beta
    X_mean = []

    for _ in range(30) :

        pop.new_deal()

        for _ in range(100)   :
            pop.evolve(10)
            X_mean += [ mean( pop.get_state_vector() ) ]

    hist( X_mean, bins = linspace(-1,1,50), alpha = 0.5, label = str( beta ), normed = True )
    plot( mean(X_mean), 0, 'o', color = 'k')
    print(mean(X_mean), r'+/-', std(X_mean)/sqrt(len(X_mean)) )

legend()
show()
