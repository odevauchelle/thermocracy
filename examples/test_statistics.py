from pylab import *
from matplotlib.tri import Triangulation

import sys
sys.path.append('./../')

from thermocracy import population, plot_connectivity, statistics, triangulation_to_connectivity, neighbors_opinion


x = rand(10)
y = rand(len(x))

Th = Triangulation(x,y)

pop = population( connectivity = triangulation_to_connectivity(Th), epsilon = 1, beta = 2, state = None )

yes_neigh, no_neigh = pop.get_neighbors_opinion()

figure()
triplot(Th)

X = pop.get_state_vector()

for i in range(len(x)) :
    if X[i] > 0 :
        plot( x[i], y[i], 'og', ms = 15, alpha = .5 )
    else :
        plot( x[i], y[i], 'or', ms = 15, alpha = .5 )
    text( x[i], y[i], str( pop.get_number_of_like_minded_neighbors()[i] ), color = 'k', ha = 'center', va = 'center')
    # text( x[i], y[i], str( (yes_neigh[i], no_neigh[i]) ), color = 'k', ha = 'center', va = 'center')

show()
