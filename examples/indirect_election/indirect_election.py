from pylab import *
from matplotlib.tri import Triangulation
from random import choice

import sys
sys.path.append('./../../')

import thermocracy as thm

x = linspace(0,1,50)
y = x
x, y = meshgrid(x,y)
x, y = x.flatten(), y.flatten()

Th = Triangulation( x, y )

N = len( Th.x )

connectivity = thm.triangles_to_connectivity( Th.triangles )

connectivity_array = connectivity.toarray()
connectivity_array += connectivity_array.T

# print(connectivity_array)
# N = len(x)

# H = thm.Hamiltonian( terms = [ thm.neighbors_influence ], coeffs = [ 1 ] )

# pop = thm.population( connectivity = connectivity, H = H, beta = None, state = None )

figure()
ax_nodes = gca()

ax_nodes.axis('equal')
ax_nodes.axis('off')
ax_nodes.triplot(Th, lw = .5, color = 'k', alpha = .2)

states = []
i = 0
free_nodes = arange( N ).tolist()
state_max_size = 100

while free_nodes != [] :

    state = [ free_nodes.pop( choice( arange( len( free_nodes ) ) ) ) ] # start a new state

    i_s = 0

    while i_s < len( state ) and len(state ) <= state_max_size :

        for i in connectivity_array[ :, state[ i_s ] ].nonzero()[0] :

            try : # collect neighbors
                state += [ free_nodes.pop( free_nodes.index( i ) ) ]
            except :
                pass
        
        i_s += 1

    states += [state]


#######################
#
# map
#
#######################


for state in states :
    ax_nodes.plot( Th.x[ state ], Th.y[ state ], '.' )

# ax_nodes.plot(x,y,'o')
# p, = ax_nodes.plot(x,y,'o')
# figure()
# ax_hist = gca()


# show(block = False)

# for pop.beta in array( [.1, .3, 1.] ) :

#     ax_nodes.set_title(r'$\beta=$' + str(pop.beta))
#     X_mean = []

#     for _ in range(20) :

#         pop.new_deal()
#         pop.evolve( 300*pop.N )

#         positive = pop.state > 0
#         p.set_data( x[positive], y[positive] )
#         pause(.01)
#         X_mean += [ mean( pop.state ) ]

#     ax_hist.hist( X_mean, alpha = 0.5, label = str( pop.beta ) , bins = linspace(-1,1,( len(X_mean)//2 )*2 + 1 ) )

# ax_hist.legend()
show()
