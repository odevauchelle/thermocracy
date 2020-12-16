from pylab import *
from matplotlib.tri import Triangulation

import sys
sys.path.append('./../')

import thermocracy as thm

x = linspace(0,1,10)
y = x
x, y = meshgrid(x,y)
x, y = x.flatten(), y.flatten()

Th = Triangulation(x,y)

connectivity = thm.triangles_to_connectivity( Th.triangles )

N = len(x)

H = thm.Hamiltonian( terms = [ thm.neighbors_influence ], coeffs = [ 1 ] )

pop = thm.population( connectivity = connectivity, H = H, beta = None, state = None )

figure()
ax_nodes = gca()
ax_nodes.plot(x,y,'o')
p, = ax_nodes.plot(x,y,'o')
ax_nodes.axis('equal')
ax_nodes.axis('off')
ax_nodes.triplot(Th, lw = .5, color = 'k', alpha = .2)

figure()
ax_hist = gca()


show(block = False)

for pop.beta in array( [.1, .3, 1.] ) :

    ax_nodes.set_title(r'$\beta=$' + str(pop.beta))
    X_mean = []

    for _ in range(20) :

        pop.new_deal()
        pop.evolve( 300*pop.N )

        positive = pop.state > 0
        p.set_data( x[positive], y[positive] )
        pause(.01)
        X_mean += [ mean( pop.state ) ]

    ax_hist.hist( X_mean, alpha = 0.5, label = str( pop.beta ) , bins = linspace(-1,1,( len(X_mean)//2 )*2 + 1 ) )

ax_hist.legend()
show()
