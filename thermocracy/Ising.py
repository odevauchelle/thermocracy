from scipy.sparse import csr_matrix
from random import getrandbits, randint
from scipy import array, rand, mean, exp, vectorize

if __name__ == '__main__' :
    from graphics import plot_connectivity
    from statistics import *
else :
    from .graphics import plot_connectivity
    from .statistics import *

##############################
#
# useful functions
#
##############################

def int_to_state_vector( state_int, N ):
    state_bin = bin(state_int)[2:]
    return array( [-1]*( N - len( state_bin ) ) +  [ 2*int(bit) - 1 for bit in state_bin ] )

def state_vector_to_int( state_vect ):
    return int( '0b' + ''.join( [ str( ( spin + 1 )//2 ) for spin in state_vect ] ), 2 )

def default_acceptance_probability( dE, beta ) :
    '''
    Metropolis-Hasting acceptance probability.
    a = min( [ 1, exp( -dE*beta ) ] )
    '''
    return min( [ 1, exp( -dE*beta ) ] )

##########################
#
# Population class
#
##########################

class population :

    def __init__( self, connectivity, H, beta, state = None, acceptance_probability = None, **kwargs ) :
        '''
        state is an integer
        '''

        self.connectivity = csr_matrix( connectivity )
        self.H = H
        self.beta = beta
        self.N = self.connectivity.shape[0]

        if acceptance_probability is None :
            self.acceptance_probability = default_acceptance_probability

        if state is None :
            self.new_deal()

        else :
            self.state = state

    def plot_connectivity( self, *args, **kwargs ) :
        plot_connectivity( self.connectivity, *args, **kwargs   )

    def new_deal( self, opinion = None ) :

        if opinion is None :
            self.state = getrandbits( self.N )
        else :
            state = ( rand( self.N ) < opinion )*2 - 1
            self.set_state( state )

    def get_state_vector( self ) :
        return int_to_state_vector( self.state, self.N )

    def set_state( self, state_vector ) :
        self.state = state_vector_to_int( state_vector )

    def flip( self, i ) :
        state_vector = self.get_state_vector()
        state_vector[i] *= -1
        self.set_state(state_vector)

    def get_opinion( self, state = None ) :

        if state is None :
            state = self.state

        return len( bin( state )[2:].replace('0','') )

    def get_E_terms( self, X = None ) :

        if X is None :
            X = self.get_state_vector()

        return self.H.get_contributions( X = X, connectivity = self.connectivity )

    def get_E( self, X = None ) :

        if X is None :
            X = self.get_state_vector()

        return self.H.get_energy( X = X, connectivity = self.connectivity )

    def evolve( self, step_number = 1 ) :

        for _ in range( step_number ) :

            # pick a node
            i = randint( 0, self.N - 1 )

            # calculate energy jump associated to flipping this node

            dE = self.H.get_energy( X = self.get_state_vector(), i = i, connectivity = self.connectivity )

            if rand() < self.acceptance_probability( dE, self.beta ) :
                # keep
                self.flip(i)

    def get_neighbors_opinion( self ) :
        return neighbors_opinion( self.connectivity, self.get_state_vector() )

    def get_number_of_neighbors( self ) :
        return number_of_neighbors( self.connectivity, self.get_state_vector() )

    def get_number_of_like_minded_neighbors( self ) :
        return number_of_like_minded_neighbors( self.connectivity, self.get_state_vector() )

if __name__ == '__main__' :

    # print(help(randint))

    from pylab import *
    from Hamiltonian import *

    C = ( rand(*[5]*2) > .5 )*1

    H = Hamiltonian( terms = [neighbors_influence, polls_influence ], coeffs = [ 1, 1 ] )

    pop = population( connectivity = C, H = H, beta = 1, state = None )

    X = pop.get_state_vector()

    for term in [neighbors_influence, polls_influence ] :
        print( term.get_contribution( X, 1, connectivity = C ) )


    print(pop.connectivity.toarray())
    print(pop.N)

    flip_i = 2
    print(pop.get_state_vector(), pop.get_E())
    pop.flip(flip_i)
    print(pop.get_state_vector(), pop.get_E())


    dE = linspace(-1/3,1,101)*6
    p = [ pop.acceptance_probability( dE_single, beta = 1 ) for dE_single in dE ]
    plot(dE, p )
    title( pop.acceptance_probability.__doc__.split('\n')[1] )
    fill_between( dE, p, alpha = .1 )
    show()
