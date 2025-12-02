#
#   plot_guess                       
#                               
#   Script to plot an initial flowfield guess created using the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_guess.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the guess data
    filename = 'out_guess_' + sys.argv[-1] + '.bin'

    # Read the case from file
    g = read_case(filename)

    # Open figure window and open four subplots
    fig,ax = plt.subplots(2,2,sharex=True,sharey=True,figsize=[14.4,7.2]); 
    fig.tight_layout()

    # Set subplot aspect ratios as equal and remove axes labels
    ax = ax.flatten()
    for a in ax:
        a.set_aspect('equal',adjustable='box'); a.axis('off')

    # Plot the primary flow variables to show the guess
    fieldnames = ['ro','roe','rovx','rovy']
    
    for n,name in enumerate(fieldnames):
        mina = []
        maxa = []
        for m in range(g['nn']):
            b = cut_block(g,m)
            mina.append( np.min(b[name]) )
            maxa.append( np.max(b[name]) )
        print(name)
        print(mina)
        print(maxa)

        minv = np.min(mina)
        maxv = np.max(maxa)

        for m in range(g['nn']):
            # Plot filled contour levels
            b = cut_block(g,m)
            hc = ax[n].pcolormesh(b['x'],b['y'],b[name],shading='gouraud', vmin = minv, vmax = maxv)

            # Draw the walls of the block
            plot_wall(ax[n],b)

  	# Add colorbar with variable name
        colorbar(hc,name)
        #scaling the limits of the graph
        xmin = min(np.min(g['x'][m]) for m in range(g['nn']))
        xmax = max(np.max(g['x'][m]) for m in range(g['nn']))
        ymin = min(np.min(g['y'][m]) for m in range(g['nn']))
        ymax = max(np.max(g['y'][m]) for m in range(g['nn']))

        ax[n].set_xlim(xmin, xmax)
        ax[n].set_ylim(ymin, ymax)
       

    # Show all the plots
    plt.savefig("myplot2.png")



    
main()


