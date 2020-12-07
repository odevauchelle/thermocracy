from pylab import *
from matplotlib.tri import Triangulation

import sys
sys.path.append('./../')

import thermocracy as thm

N = 100


H = thm.Hamiltonian( terms = [ thm.neighbors_influence, thm.polls_influence], coeffs = [ nan, nan ] )

# theta = linspace(0, 2*pi, N+1)[:-1]
# x = cos(theta)
# y = sin(theta)
#
# edges = [ [i,i+1] for i in range(len(x)-1) ]
# edges += [ [ len(x) - 1, 0 ] ]
#
x = linspace( 0,1, 10 )
y = x
x,y = meshgrid(x,y)
x = x.flatten()
y = y.flatten()

Th = Triangulation(x,y)
edges = Th.edges



connectivity = thm.edges_to_connectivity( edges, size = len(x) )

figure()
ax_nodes = gca()
ax_nodes.axis('equal')

for edge in edges :
    ax_nodes.plot( x[edge], y[edge], color = 'k', alpha = .2 )

figure()
ax_hist = gca()


pop = thm.population( connectivity = connectivity, H = H, beta = None, state = None )

pop.beta = 100
evolve_step = 10

show(block = False)

for epsilon in [ 1, 5 ] :

    pop.H.set_coeffs( [ 1, epsilon ] )

    X_mean = []

    for _ in range(30) :

        pop.new_deal()

        for _ in range( int( 5*N/evolve_step ) )  :
            pop.evolve( evolve_step )
            X_mean += [ mean( pop.get_state_vector() ) ]

    ax_hist.hist( X_mean, bins = linspace( -1, 1, int( len( X_mean )/30 ) ), alpha = 0.5, label = str( pop.H.coeffs['polls'] ), normed = True )
    yes = pop.get_state_vector() > 0
    ax_nodes.plot( x[yes], y[yes], '.', color = 'tab:blue' )
    ax_nodes.plot( x[~yes], y[~yes], '.', color = 'tab:orange' )

    print(mean(X_mean), r'+/-', std(X_mean)/sqrt(len(X_mean)) )

    pause(0.01)
    input('?')

legend()
show()
