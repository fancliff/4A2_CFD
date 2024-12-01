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

    # For the sod shock tube we want to plot velocity ratio relative to v_ref
    # We also want to plot the density ratio relative to ro_ref
    # Both of these quantities can be plotted at j_mid = (nj+1) // 2
    # Could also plot the quantities along either both or either of
    #   the uper and lower walls.

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
    g = calc_secondary(av,g)    


    g['ro_ratio'] = g['ro']/1.0
    g['v_ratio'] = g['v']/2.846049894


    # Specify the parameters to plot
    fieldnames = ['ro_ratio','v_ratio']; 
    colnames = ['Density Ratio','Velocity ratio']

    j_mid = (g['nj'])//2
    # j_pos = 0
    j_pos = j_mid

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):

        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box')
 
        # Plot the quantity against x for t =0.2s and 0.0s
        ax.plot(g['x'][:,j_pos],g[name][:,j_pos],label='t=0.2s',color='blue')
        if name == 'ro_ratio':
            t0 = np.ones(g['ni'])
            t0[(g['ni']//2):g['ni']] *= 0.125
        if name == 'v_ratio':
            t0 = np.zeros(g['ni'])
        ax.plot(g['x'][:,j_pos],t0,label='t=0.0s',color='red')

        # Add labels and title (fortran indexing for j_pos printout)
        ax.set_title(f'{colnames[n]} against x at j = {j_pos+1}', fontsize = 14)
        ax.set_xlabel('x', fontsize = 12)
        ax.set_ylabel(colnames[n], fontsize = 12)
        ax.legend(loc = 'best', fontsize = 12)

        ax.set_ylim(-0.05,1.05)
        ax.set_xlim(-0.05,1.05)

    # Show all the plots
    plt.show()

    
main()


