from pylab import *

import sys
sys.path.append('./../')

import thermocracy as thm

nb_of_nodes = 50
connectivity = ( rand( *[nb_of_nodes]*2 ) > .9 )*1

H = thm.Hamiltonian( terms = [ thm.neighbors_influence, thm.polls_influence ], coeffs = [ 1, 1 ] )
pop = thm.population( connectivity = connectivity, H = H, beta = .2, state = None )

X_mean = [ mean( pop.state ) ]

for _ in range(300)   :
    pop.evolve(10)
    X_mean += [ mean( pop.state ) ]

theta =  linspace( 0, 2*pi, nb_of_nodes + 1 )[:-1]
x, y = cos(theta), sin(theta)



figure()
thm.plot_connectivity( connectivity = pop.connectivity, x = x, y = y, color = 'LightGrey' )
plus = pop.state > 0
plot( x[plus], y[plus],'o', color = 'tab:blue')
plot( x[~plus], y[~plus],'o', color = 'tab:orange')

axis('equal')
axis('off')
# savefig('Connectivity.svg', bbox_inches = 'tight')

figure()
plot( X_mean )
xlabel('time')
ylabel('Magnetization')
# savefig('Magnetization.svg', bbox_inches = 'tight')

show()
