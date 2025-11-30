#
#   plot_coord
#                               
#   Script to plot a mesh created using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_coord.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filename to read the grid data
    filename = 'out_coord_' + sys.argv[-1] + '.bin'

    # Read the case from file
    g = read_case(filename)

    # Open figure window and set the axes to be equal
    plt.figure(figsize=[9.6,7.2]); ax = plt.axes(); cols = gen_cols();
    ax.set_xlabel('x / m'); ax.set_ylabel('y / m');
    ax.set_aspect('equal',adjustable='box'); ax.tick_params(direction='in')

    for m in range(g['nn']):

        # Plot the mesh coordinates to show the cells
        ax.plot(g['x'][m],g['y'][m],color=cols[0,:],linewidth=0.5)  #change the colours per block
        ax.plot(np.transpose(g['x'][m]),np.transpose(g['y'][m]),color=cols[0,:],
            linewidth=0.5)

    # Draw the boundary of the block
        plot_bound(ax,cut_block(g,m))

    #scaling the limits of the graph
    xmin = min(np.min(g['x'][m]) for m in range(g['nn']))
    xmax = max(np.max(g['x'][m]) for m in range(g['nn']))
    ymin = min(np.min(g['y'][m]) for m in range(g['nn']))
    ymax = max(np.max(g['y'][m]) for m in range(g['nn']))

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    # Show all the plots
    plt.savefig("myplot.png")

    
main()


