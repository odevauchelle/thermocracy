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
    Metropolis-Hasting acceptance probability
    a = min( [ 1, exp( -dE*beta ) ] )
    '''
    return min( [ 1, exp( -dE*beta ) ] )

##########################
#
# Population class
#
##########################

class population :

    def __init__( self, connectivity, H = None, beta = None, state = None, acceptance_probability = None, **kwargs ) :
        '''
        state is a list of -1 and 1
        '''

        self.connectivity = csr_matrix( connectivity )
        self.H = H
        self.beta = beta
        self.N = self.connectivity.shape[0]
        self.E = None

        if acceptance_probability is None :
            self.acceptance_probability = default_acceptance_probability

        if state is None :
            self.new_deal()

        else :
            try :
                self.state = state_vector_to_int( state )
            except :
                self.state = state

    def plot_connectivity( self, *args, **kwargs ) :
        plot_connectivity( self.connectivity, *args, **kwargs   )

    def new_deal( self, magnetization = None ) :

        if magnetization is None :
            magnetization = 0.

        self.set_state( ( rand( self.N ) < ( 1 + magnetization )/2 )*2 - 1 )

    def get_state( self ) :
        return self.state

    def get_integer_state( self ) :
        return state_vector_to_int( self.state )

    def set_state( self, state ) :

        try :
            state[0]
            self.state = state

        except :
            self.state = int_to_state_vector( state, self.N )

        self.E = None

    def flip( self, i ) :
        self.state[i] *= -1
        self.E = None

    def get_opinion( self, state = None ) :

        if state is None :
            state = self.state

        return mean( state )

    def get_E_terms( self, X = None ) :

        if X is None :
            X = self.state

        return self.H.get_contributions( X = X, connectivity = self.connectivity )

    def get_E( self, X = None ) :

        if X is None :
            X = self.state

        self.E = self.H.get_energy( X = X, connectivity = self.connectivity )
        return self.E

    def evolve( self, step_number = 1 ) :

        for _ in range( step_number ) :

            if self.E is None :
                self.get_E()

            E_old = self.E

            # pick a node
            i = randint( 0, self.N - 1 )

            # flip it

            self.flip(i)

            # calculate new energy

            self.get_E()
            dE = self.E - E_old

            if rand() < self.acceptance_probability( dE, self.beta ) :
                # accept
                pass
            else :
                # reject
                self.flip(i)
                self.E = E_old


    def get_neighbors_opinion( self ) :
        return neighbors_opinion( self.connectivity, self.get_state_vector() )

    def get_number_of_neighbors( self ) :
        return number_of_neighbors( self.connectivity, self.get_state_vector() )

    def get_number_of_like_minded_neighbors( self ) :
        return number_of_like_minded_neighbors( self.connectivity, self.get_state_vector() )

if __name__ == '__main__' :

    from pylab import *

    print('Test state vector to int:')

    sqrt_N = 100
    N = sqrt_N**2
    print('N =',N)
    X = 2*( rand(N) > .5 ) - 1

    imshow(X.reshape([sqrt_N]*2))

    X_c = int_to_state_vector( int( str( state_vector_to_int( X ) ) ), N )
    print( sum( X-X_c )==0 )

    print('------------------------')



    from Hamiltonian import *

    C = ( rand(*[5]*2) > .5 )*1

    H = Hamiltonian( terms = [ neighbors_influence, polls_influence ], coeffs = [ 1, 1 ] )

    pop = population( connectivity = C, H = H, beta = 1, state = None )

    X = pop.state
    print(X)
    print(pop.N)

    for term in [neighbors_influence, polls_influence ] :
        print( term.get_contribution( X, 1, connectivity = C ) )


    print(pop.connectivity.toarray())
    print(pop.N)

    flip_i = 2
    print(pop.state, pop.get_E())
    pop.flip(flip_i)
    print(pop.state, pop.get_E())


    dE = linspace(-1/3,1,101)*6
    p = [ pop.acceptance_probability( dE_single, beta = 1 ) for dE_single in dE ]
    figure()
    plot(dE, p )
    title( pop.acceptance_probability.__doc__.split('\n')[1] )
    fill_between( dE, p, alpha = .1 )
    show()
