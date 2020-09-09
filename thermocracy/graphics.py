
from pylab import *
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

def merge( default, new ) :

    if new is None :
        return default

    else :
        intermediate = default.copy()
        intermediate.update(new)
        return intermediate

default_frame_style = dict( facecolor = 'none', clip_on = False )
default_matrix_style = dict( complementary_color = 'w' )
default_text_style = dict( color = 'inverted', va = 'center', ha = 'center' )

def plot_connectivity( C, ax = None, frame = True, indices = True, frame_style = None, matrix_style = None, text_style = None ) :

    if ax is None :
        ax = gca()

    matrix_style = merge( merge( default_matrix_style, dict( color = next(ax._get_lines.prop_cycler)['color']  ) ), matrix_style )
    frame_style = merge( merge( default_frame_style, dict( edgecolor = matrix_style['color'] ) ), frame_style )
    text_style = merge( default_text_style, text_style )

    N = C.shape[0]

    try :
        C = C.toarray()
    except :
        pass

    cmap = mpl.colors.ListedColormap([matrix_style['complementary_color'], matrix_style['color']])
    bounds = [-6,.5,6]
    norm = mpl.colors.BoundaryNorm( bounds, cmap.N )

    ax.imshow( C, norm = norm, cmap = cmap )

    if frame :
        ax.add_collection( PatchCollection( [ Rectangle( [-.5, -.5 ], N, N ) ], **frame_style ) )

    if text :

        invert = text_style['color'] == 'inverted'

        for i in range( N ) :

            ij = i, 0

            zero = True

            for ij in ij, ij[::-1] :

                if invert :

                    if C[ ij[1], ij[0] ] > 0.5 :
                        text_style['color'] = matrix_style['complementary_color']
                    else :
                        text_style['color'] = matrix_style['color']

                if zero or sum(ij) > 0 :
                    if sum(ij) == 0 :
                        zero = False

                    text( *ij, str(i), **text_style )

    ax.axis("equal")
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])

if __name__ == '__main__' :

    C = ( rand(8,8) > .5 )*2 - 1

    figure()

    plot_connectivity(C)

    show()
