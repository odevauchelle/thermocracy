from scipy.sparse import csr_matrix
from random import getrandbits, randint
from scipy import array, rand, mean, exp

def int_to_state_vector( state_int, N ):
    state_bin = bin(state_int)[2:]
    return array( [-1]*( N - len( state_bin ) ) +  [ 2*int(bit) - 1 for bit in state_bin ] )

def state_vector_to_int( state_vect ):
    return int( '0b' + ''.join( [ str( ( spin + 1 )//2 ) for spin in state_vect ] ), 2 )

def H( X, M, epsilon, N ):
    return -( M.dot( X ) ).dot( X )/N + epsilon*mean( X )**2

def dH( X, M, epsilon, N, i ) :
    E_old = H( X, M, epsilon, N )
    X[i] *= -1
    return H( X, M, epsilon, N ) - E_old

def probability( dE, beta ) :
    p_Boltzmann = exp( -dE*beta )
    return p_Boltzmann/( 1. + p_Boltzmann )

class population :

    def __init__( self, connectivity, epsilon, beta, state = None ) :
        '''
        state is an integer
        '''

        self.connectivity = csr_matrix( connectivity )
        self.epsilon = epsilon
        self.beta = beta
        self.N = self.connectivity.shape[0]

        if state is None :
            self.new_deal()

        else :
            self.state = state

    def new_deal( self ) :
        self.state = getrandbits( self.N )

    def get_state_vector( self ) :
        return int_to_state_vector( self.state, self.N )

    def set_state( self, state_vector ) :
        self.state = state_vector_to_int( state_vector )

    def flip( self, i ) :
        bin_state = bin( self.state )[2:]
        bin_state = '0'*( self.N - len(bin_state) ) + bin_state
        self.state = int( '0b' + bin_state[:i] + str( int( not( int( bin_state[i] ) ) ) ) + bin_state[i+1:], 2 )

    def get_E( self, X = None ) :

        if X is None :
            return H( self.get_state_vector(), self.connectivity, self.epsilon, self.N )

        else :
            return H( X, self.connectivity, self.epsilon, self.N )

    def evolve( self, step_number = 1 ) :

        # initialize

        X = self.get_state_vector().copy()
        E = self.get_E()

        for _ in range( step_number ) :

            # pick a node
            n = randint(0,self.N-1)

            # flip it
            X[n] *= -1

            # get new energy:
            E_new = self.get_E( X )

            # keep or drop according to Boltzmann

            if rand() < probability( E_new - E, self.beta ) :
                # keep
                E = E_new

            else :
                # drop
                X[n] *= -1

        # reccords
        self.set_state( X )



if __name__ == '__main__' :

    # print(help(randint))

    from pylab import *

    C = ( rand(*[5]*2) > .5 )*1

    pop = population( connectivity = C, epsilon = 1, beta = 1, state = None )

    print(pop.connectivity.toarray())
    print(pop.N)

    flip_i = 2
    print(pop.get_state_vector(), pop.get_E())
    pop.flip(flip_i)
    print(pop.get_state_vector(), pop.get_E())
