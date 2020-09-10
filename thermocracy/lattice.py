
from scipy.sparse import csr_matrix
from numpy import array

def triangle_to_edges( triangle ) :
    return [ tuple( sorted( [ triangle[ij[0]], triangle[ij[1]] ] ) ) for ij in [ [0,1], [1,2], [2,0] ] ]

def triangles_to_edges( triangles ) :

    edges = []

    for triangle in triangles :
        edges += triangle_to_edges( triangle )

    return list( set( edges ) )

def triangles_to_node_indices( triangles ) :

    return list( set( triangles.flatten() ) )

def triangles_to_connectivity( triangles, size = None ) :

    edges = triangles_to_edges( triangles )

    if size is None :
        size = len( triangles_to_node_indices( triangles ) )

    return csr_matrix( ( [1]*len( edges ), array( edges ).T ), shape = [ size ]*2 ).transpose()

def triangulation_to_connectivity( Th ) :
    return csr_matrix( ( [1]*len( Th.edges ), Th.edges.T ), shape = [ len( Th.x ) ]*2 )


if __name__ == '__main__' :

    from pylab import *
    from matplotlib.tri import Triangulation

    x = rand(6)
    y = rand(len(x))
    Th = Triangulation( x, y )

    figure('mesh')
    triplot(Th)
    axis('equal')

    edges = triangles_to_edges( Th.triangles )

    M1 = triangulation_to_connectivity( Th )
    M2 = triangles_to_connectivity( Th.triangles )

    print(M1.toarray())
    print(M2.toarray())

    show()
