
from scipy.sparse import csr_matrix
from scipy.sparse import triu as sparse_triu
from numpy import array, nan
from matplotlib.pyplot import gca
from networkx import adjacency_matrix, kamada_kawai_layout

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

def connectivity_to_edges( connectivity ) :

    row, col = connectivity.nonzero()

    return [ [ row[i], col[i] ] for i in range( len(row) ) ]

def triangles_to_connectivity( triangles, size = None ) :

    edges = triangles_to_edges( triangles )

    if size is None :
        size = len( triangles_to_node_indices( triangles ) )

    return edges_to_connectivity( edges, size )

def normalize_connectivity( connectivity, axis = 0 ) :

    connectivity += connectivity.T # in case input is not symmetric

    connectivity = connectivity.multiply( 1./connectivity.sum( axis = axis ) )

    connectivity += connectivity.T

    connectivity /= connectivity.sum()/connectivity.count_nonzero()

    return connectivity


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

def plot_edges( x, y, edges, ax = None, **kwargs ) :

    x_p = []
    y_p = []

    for edge in edges :
        i, j = edge
        x_p += [ x[i], x[j], nan ]
        y_p += [ y[i], y[j], nan ]

    if ax is None :
        ax = gca()

    return ax.plot( x_p, y_p, **kwargs )



def plot_connectivity( show_link_strength = True, **kwargs ):
    '''
    Convenience function. Plots connectivity in the form of triangulation, connectivity matrix or edges.
    '''

    try :
        ax = kwargs.pop('ax')
    except :
        ax = gca()

    try :
        triangulation = kwargs.pop('triangulation')
        return ax.triplot( triangulation, **kwargs )

    except :

        x = kwargs.pop('x')
        y = kwargs.pop('y')

        try :
            triangles = kwargs.pop('triangles')
            return ax.triplot( x, y, triangles, **kwargs )

        except :

            try :
                edges = kwargs.pop('edges')
                return plot_edges( x, y, edges, ax = ax, **kwargs )

            except :

                C = kwargs.pop('connectivity')
                C_max = C.max()

                if show_link_strength :

                    x = array(x)
                    y = array(y)

                    try :
                        base_line_width = kwargs.pop( 'linewidth' )
                    except :
                        base_line_width = rcParams['lines.linewidth']

                    try :
                        base_alpha = kwargs.pop( 'alpha' )
                    except :
                        base_alpha = 1

                    try :
                        color = kwargs.pop('color')
                    except :
                        color = 'tab:blue'

                    for i in range( C.get_shape()[0] ) :
                        for j in range(i) : # assumes matrix is symmetric
                            if C[i,j] != 0 :
                                ax.plot( x[[i,j]], y[[i,j]], color = color, lw = C[i,j]*base_line_width, alpha = C[i,j]*base_alpha/C_max, **kwargs  )

                    return None

                else :
                    edges = connectivity_to_edges( connectivity )
                    return plot_edges( x, y, edges, ax = ax, **kwargs )





def extract_mesh_data( p_mesh ) :

    for mesh_keys in [ ['triangulation'], ['x','y','triangles'], [ 'x', 'y', 'edges' ], ['x','y', 'connectivity'] ] :

        try :
            mesh_data = {}

            for key in mesh_keys :
                mesh_data[ key ] = p_mesh[ key ]

            return mesh_data

        except :
            pass


def graph_to_edges( graph, layout = None ) :

    '''
    Convert networkx graph into nodes and edges.
    '''

    if layout is None :
        layout = kamada_kawai_layout

    connectivity = sparse_triu( adjacency_matrix( graph ) )
    pos = layout(graph)

    N = connectivity.shape[0]

    x = [0]*N
    y = [0]*N

    for i in range( N ) :
        x[i], y[i] = pos[i].tolist()

    edges = array( connectivity_to_edges( connectivity ) ).tolist()

    return x, y, edges

##############################
#
# TRY IT OUT
#
##############################

if __name__ == '__main__' :

    from pylab import *
    from matplotlib.tri import Triangulation

    x = rand(15)
    y = rand(len(x))
    Th = Triangulation( x, y )

    edges = triangles_to_edges( Th.triangles )

    M1 = triangulation_to_connectivity( Th )
    M2 = triangles_to_connectivity( Th.triangles )

    print('Test normalize_connectivity')
    M = M1
    print('Original')
    print( M.toarray() )

    M = normalize_connectivity( M )
    print('Normalized')
    print( M.toarray() )
    print('Mean non-zero coefficient')
    print( mean( M.toarray()[ M.toarray() != 0 ] ) )


    print('---------------------------')



    figure('mesh')
    # plot_edges(x,y,edges)
    plot_connectivity( x = x, y = y, connectivity = M1 )
    plot_connectivity( triangulation = Th, color = 'r', linestyle = '--', alpha = .3 )


    axis('equal')


    # print(M1)
    M3 = connectivity_to_dict(M2)
    # print(type(M3['col'][0]))
    # print(dict_to_connectivity(connectivity_to_dict(M2, with_data = False )))

    figure()

    from networkx import barabasi_albert_graph

    graph = barabasi_albert_graph( n = 70, m = 2)

    x, y, edges = graph_to_edges( graph )
    C = normalize_connectivity( edges_to_connectivity( edges, size = len(x) ) )
    plot_connectivity( x = x, y = y, connectivity = C )
    axis('equal')

    print('Mean non-zero coefficient')
    print( mean( C.toarray()[ C.toarray() != 0 ] ) )

    show()
