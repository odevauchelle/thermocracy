from pylab import *
from matplotlib.tri import Triangulation
import json

import sys
sys.path.append('./../')

import thermocracy as thm


H = thm.Hamiltonian( terms = [ thm.neighbors_influence, thm.polls_influence ], coeffs = [ 1, 1 ] )

N = 4
theta = linspace(0, 2*pi, N+1)[:-1]
x = cos(theta)
y = sin(theta)
edges = [ [i,i+1] for i in range(len(x)-1) ]
edges += [ [ len(x) - 1, 0 ] ]
connectivity = thm.edges_to_connectivity( edges, size = len(x) )

figure()
ax_nodes = gca()
ax_nodes.axis('equal')
ax_nodes.axis('off')


for edge in edges :
    ax_nodes.plot( x[edge], y[edge], color = 'k', alpha = .2 )

state = array([ 1, -1, 1, -1 ])

pop = thm.population( connectivity, H, 1  )
pop.set_state( state )

print('terms',pop.get_E_terms())
print('coefs',pop.H.coeffs)
print('E',pop.get_E())
i = 3
print('dE',i,pop.H.get_energy( X = pop.get_state_vector(), i = i, connectivity = pop.connectivity))

pop.flip(i)
print('terms',pop.get_E_terms())
print('coefs',pop.H.coeffs)
print('E',pop.get_E())


yes = pop.get_state_vector() > 0
p_yes, = ax_nodes.plot( x[yes], y[yes], 'o', color = 'tab:blue' )
p_no, = ax_nodes.plot( x[~yes], y[~yes], 'o', color = 'tab:orange' )

for i,_ in enumerate(x) :
    ax_nodes.text( x[i], y[i], str(i) )

show()
