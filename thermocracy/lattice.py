
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

    return list( set( array( triangles ).flatten() ) )

def edges_to_connectivity( edges, size ) :
    return csr_matrix( ( [1]*len( edges ), array( edges ).T ), shape = [ size ]*2 ).transpose()

def triangles_to_connectivity( triangles, size = None ) :

    edges = triangles_to_edges( triangles )

    if size is None :
        size = len( triangles_to_node_indices( triangles ) )

    return edges_to_connectivity( edges, size )

def triangulation_to_connectivity( Th ) :
    return csr_matrix( ( [1]*len( Th.edges ), Th.edges.T ), shape = [ len( Th.x ) ]*2 )

def connectivity_to_dict( connectivity, with_data = True ) :

    connectivity = connectivity.tocoo()

    dic = dict(
        shape = array( connectivity.shape ).tolist(),
        row = connectivity.row.tolist(),
        col = connectivity.col.tolist()
        )

    if with_data :
        dic['data'] = connectivity.data.tolist()

    return dic

def dict_to_connectivity( dic ) :
    try :
        data = dic['data']
    except :
        data = [1]*len(dic['row'])

    return csr_matrix( ( data, (dic['row'], dic['col'])), shape = dic['shape']  )

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

    print(M1)
    M3 = connectivity_to_dict(M2)
    print(type(M3['col'][0]))
    print(dict_to_connectivity(connectivity_to_dict(M2, with_data = False )))

    show()
