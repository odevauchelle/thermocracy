from pylab import *
from matplotlib.tri import Triangulation
from random import choice
from scipy.sparse import csr_matrix

import sys
sys.path.append('./../../')

import thermocracy as thm

x = linspace(0,1,50)
y = x
x, y = meshgrid(x,y)
x, y = x.flatten(), y.flatten()

Th = Triangulation( x, y )

N = len( Th.x ) # electorate size

connectivity = thm.triangles_to_connectivity( Th.triangles )

connectivity_array = connectivity.toarray()
connectivity_array += connectivity_array.T

####################
# 
# Create states
#
####################

states = []
state_size = []
i = 0
free_nodes = arange( N ).tolist()
state_max_size = 200

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
    state_size += [len(state)]

state_size = array(state_size)

N_s = len( state ) 

#####################
#
# Get state borders
#
#####################

state_index = x*0
i = 0

for i, state in enumerate(states ) :
    state_index[state] = i

boundary = []

for triangle in Th.triangles :
    if len( set( state_index[ triangle ] ) ) != 1 :
        boundary += [ [ mean( Th.x[triangle] ), mean( Th.y[triangle]) ] ]

boundary = array(boundary)

#####################
#
# figures
#
#####################


fig, ( ax_nodes, ax_hist ) = subplots( ncols = 2, figsize = array([ 1,1/2 ])*15 )

for ax in ( ax_nodes, ) :
    ax.axis('equal')
    ax.axis('off')
    # ax.triplot(Th, lw = .5, color = 'k', alpha = .2)
    ax.plot( *boundary.T, '.k' )

show(block = False)
pause(0.01)

#######################
#
# map
#
#######################

# for state in states :
#     ax_states.plot( Th.x[ state ], Th.y[ state ], '.' )

# pause(0.01)

#######################
#
# Hamiltonian
#
#######################

state_result_matrix = []

for state in states :
    state_line = array( [0]*N )
    state_line[ state ] = 1
    state_result_matrix += [  state_line ]

state_result_matrix = csr_matrix( array( state_result_matrix ) )

def get_state_result( electorate_state ) :
    return state_result_matrix@electorate_state

state_result_influence = thm.HTerm(
    name = 'state_result',
    function_E = lambda X, connectivity, **kwargs: mean( sign( get_state_result( X ) ) )*sum( X )
    )   

epsilon = 0.3
H = thm.Hamiltonian( terms = [ thm.neighbors_influence, state_result_influence ], coeffs = [ 1, epsilon ] )

#######################
#
# Run simulation
#
#######################

pop = thm.population( connectivity = connectivity, H = H, beta = None, state = None )

nodes_kwargs = dict( marker = 'o', linestyle = 'none', ms = 10, alpha = .3 )

p_plus, = ax_nodes.plot( x, y, **nodes_kwargs)
p_minus, = ax_nodes.plot( x, y, **nodes_kwargs )

show(block = False)

state_result = []
country_result = []

pop.beta = 1
nb_runs = 10
bins = linspace(-1,1,10)

for i in range( nb_runs ) :

    print( round(i/nb_runs,2)*100, '%' )

    pop.new_deal()
    pop.evolve( 300*pop.N )

    positive = pop.state > 0
    p_plus.set_data( x[positive], y[positive] )
    p_minus.set_data( x[~positive], y[~positive] )

    state_result += ( get_state_result( pop.state )/state_size ).tolist()
    country_result += [ mean( sign( get_state_result( pop.state ) ) ) ]

    ax_hist.cla()
    ax_hist.hist( state_result, alpha = .5, bins = bins )
    ax_hist.set_xlabel('State electoral result')

    # ax_country.cla()
    # ax_country.hist( country_result, alpha = .5, bins = bins )
    # ax_country.set_xlabel('Country-wide result')
    # ax_country.set_xlim( [ -1, 1 ] )        

    pause(.01)

savefig('states2.pdf')

show()
