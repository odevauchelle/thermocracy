from numpy import array, shape

def opinion_to_percentage(opinion) :

    return ( opinion + 1 )*50

def neighbors_opinion( connectivity, state ) :

    yes = ( state > 0 )*1
    no = ( state < 0 )*1

    yes_neighbors = ( connectivity + connectivity.T )*yes
    no_neighbors = ( connectivity + connectivity.T )*no

    return yes_neighbors, no_neighbors

def number_of_neighbors( connectivity ) :

    return ( connectivity + connectivity.T )*array( [1]*shape(connectivity)[0] )

def number_of_like_minded_neighbors( connectivity, state ) :

    yes_neighbors, no_neighbors = neighbors_opinion( connectivity, state )

    yes = ( state > 0 )*1
    no = ( state < 0 )*1

    return yes_neighbors*yes + no_neighbors*no

def herding_ratio( connectivity, state ) :
    return number_of_like_minded_neighbors( connectivity, state )/number_of_neighbors( connectivity )

def energy_density( H, connectivity, state ) :
    return H.get_energy( connectivity = connectivity, X = state )/len(state)
