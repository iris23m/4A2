#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *

def main():

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    outname = 'out_final_' + sys.argv[-1] + '.bin'

    # Read the settings and the case from file
    av = read_settings(inname)
    g = read_case(outname)

    # When presenting results all values should be non-dimensionalised. Two
    # variables of interest might be:
    #    1. Static pressure coefficient, (p - p_ref) / (pstag_ref - p_ref)
    #    2. Mach number, v / (ga * rgas * t)**0.5
    block = []
  
    for m in range(g['nn']):
        b = cut_block(g,m)
    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
        b = calc_secondary(av,b)    

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    # INSERT
        inletcut = cut_i(b,0)
        p_ref, mass = mass_av(inletcut, 'p')
        pstag_ref, mass = mass_av(inletcut, 'pstag')
        print(b['p'].shape)
        print(p_ref)
        print(pstag_ref)
        b['cp'] = (b['p'] - p_ref) / (pstag_ref-p_ref)
        b['cpstag'] =(b['pstag'] - pstag_ref) / (pstag_ref-p_ref)
        #g['arde'] = g['dro']/p_ref
        block.append(b)
        


    # Specify the parameters to plot
    fieldnames = ['cp', 'mach', 'cpstag']; 
    colnames = ['Static pressure coefficient','Mach number', 'Stagnation pressure coefficient']

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):
        mina = []
        maxa = []
        for m in range(g['nn']):
            mina.append( np.min(block[m][name]) )
            maxa.append( np.max(block[m][name]) )
        minv = np.min(mina)
        maxv = np.max(maxa)

        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')

        for m in range(g['nn']):
        # Plot filled contour levels
            hc = ax.pcolormesh(block[m]['x'],block[m]['y'],block[m][name],shading='gouraud', vmin=minv, vmax=maxv)
            if name == 'cp':
                print(m)
                print(block[m][name])
        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        for m in range(g['nn']):
            # Add Mach = 1 contours
            if name == 'mach':
                ax.contour(block[m]['x'],block[m]['y'],block[m]['mach'],[1.0],colors='w',
                    linewidths=0.5) 

            # Draw the walls of the block
            plot_wall(ax,block[m])
        filename = f"myplot_{name}.png" 
        plt.savefig(filename)


    # Show all the plots
    plt.show()
    
    

    
main()


