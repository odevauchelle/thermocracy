
def neighbors_opinion( connectivity, state_vector ) :

    yes = ( state_vector > 0 )*1
    no = ( state_vector < 0 )*1

    yes_neighbors = ( connectivity + connectivity.T )*yes
    no_neighbors = ( connectivity + connectivity.T )*no

    return yes_neighbors, no_neighbors

def number_of_neighbors( connectivity, state_vector ) :

    return ( connectivity + connectivity.T )*abs( state_vector )

def number_of_like_minded_neighbors( connectivity, state_vector ) :

    yes_neighbors, no_neighbors = neighbors_opinion( connectivity, state_vector )

    yes = ( state_vector > 0 )*1
    no = ( state_vector < 0 )*1

    return yes_neighbors*yes + no_neighbors*no
