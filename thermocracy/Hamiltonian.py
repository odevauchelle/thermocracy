from scipy.sparse import csr_matrix
from random import getrandbits, randint
from scipy import array, rand, mean, exp

##########################
#
# Hamiltonian term class
#
##########################

class HTerm :

    def __init__( self, name, function_E, function_dE = None ) :

        self.function_E = function_E
        self.name = name

        if function_dE is None :
            def function_dE( X, i, **kwargs ) :
                old_E = self.function_E( X = X, **kwargs )
                X[i] *= -1
                return  self.function_E( X = X, **kwargs ) - old_E

        self.function_dE = function_dE

    def get_contribution( self, X, i = None, **kwargs ) :

        if i is None :
            return { self.name: self.function_E( X, **kwargs ) }

        else :
            return { self.name: self.function_dE( X, i, **kwargs ) }

##########################
#
# Hamiltonian class
#
##########################

class Hamiltonian :

    def __init__( self, terms, coeffs ) :

        try :
            terms[0]
            self.terms = terms

        except :
            self.terms = [terms]

        try :
            coeffs[ self.terms[0] ]
            self.coeffs = coeffs

        except :
            self.coeffs = {}
            for k, term in enumerate( terms ) :
                self.coeffs[ term.name ] = coeffs[k]

    def get_contributions( self,  **kwargs ) :

        contributions = {}

        for term in self.terms :
            contributions.update( term.get_contribution( **kwargs ) )

        return contributions

    def get_energy( self, **kwargs ) :

        energy = 0.

        for name, contribution in self.get_contributions( **kwargs ).items() :
            energy += self.coeffs[name]*contribution

        return energy

##########################
#
# Standard Hamiltonian terms
#
##########################

neighbors_influence = HTerm(
    name = 'neighbors',
    function_E = lambda X, connectivity, **kwargs: -( connectivity.dot( X ) ).dot( X )
    )

polls_influence = HTerm(
    name = 'polls',
    function_E = lambda X, **kwargs: mean(X)**2*len(X)
    )

bias = HTerm(
    name = 'bias',
    function_E = lambda X, **kwargs: sum(X)
)

##########################
#
# Sandbox
#
##########################

if __name__ == '__main__' :

    from pylab import *

    C = ( rand(*[15]*2) > .5 )*1

    H = Hamiltonian( terms = [ neighbors_influence, polls_influence ], coeffs = [ 1, 1 ] )

    X = ( rand(shape(C)[0]) < .5 )*2 - 1
    print(X)

    print( H.get_contributions( X = X, i = 1, connectivity = C ) )
