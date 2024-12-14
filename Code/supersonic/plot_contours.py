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

    # First complete the "calc_secondary" function within "routines.py" to
    # calculate static pressure and Mach number, and any others you want!
    g = calc_secondary(av,g)    

    # Use the "cut_i", "mass_av" AND "area_av" functions to calculate the
    # reference pressures at the inlet plane and therefore the static pressure
    # coefficient
    g_inlet = cut_i(g,0)
    pstag_ref,_ = mass_av(g_inlet,'pstag')
    p_ref,_ = area_av(g_inlet,'p')
    tstag_ref,mdot_ref = mass_av(g_inlet,'tstag')
    mdot_ref = np.sum(mdot_ref)
    #tstag_ref = av['tstag']

    g_outlet = cut_i(g,g['ni']-1)
    tstag_out,mdot_out = mass_av(g_outlet,'tstag')
    mdot_out = np.sum(mdot_out)
    tstag_out_in = tstag_out/tstag_ref
    mdot_out_in = mdot_out/mdot_ref
    
    print('')
    print(f'Reference static pressure: {p_ref:.0f}, Reference stagnation pressure: {pstag_ref:.0f}\n')
    print(f'T0exit/T0in: {tstag_out_in:.4f}')
    print(f'Mdot_exit/Mdot_in: {mdot_out_in:.4f}\n')

    g['cp'] = (g['p'] - p_ref)/(pstag_ref-p_ref)
    g['cpstag'] = (g['pstag'] - pstag_ref)/(pstag_ref-p_ref)

    # Specify the parameters to plot
    fieldnames = ['cp', 'mach','cpstag']; 
    colnames = ['Static pressure coefficient','Mach number','Stagnation pressure coefficient']

    # Plot the calculated non-dimensional parameters to show the flow solution
    for n,name in enumerate(fieldnames):

        # Open figure window
        fig = plt.figure(figsize=[9.6,7.2]); ax = plt.axes();
    
        # Set aspect ratio as equal and remove axes labels
        ax.set_aspect('equal',adjustable='box'); ax.axis('off')
 
        # Plot filled contour levels
        hc = ax.pcolormesh(g['x'],g['y'],g[name],shading='gouraud')

        # Add colorbar with variable name
        colorbar(hc,colnames[n])

        # Add Mach = 0.1 contours from mach 1 onwards
        levels = np.arange(1,3,0.1)
        if name == 'mach':
            contour = ax.contour(g['x'],g['y'],g['mach'],levels=levels,
                                 colors='w',linewidths=0.5)
            #ax.clabel(contour, inline=False, fontsize=8)
            
        # Add Pstag = 0.05 contours
        levels = np.arange(-2,2,0.02)
        if name == 'cpstag':
            contour = ax.contour(g['x'],g['y'],g['cpstag'],levels=levels,
                                 colors='w',linewidths=0.5)
            ax.clabel(contour, inline=False, fontsize=8)

        # Draw the walls of the block
        plot_wall(ax,g)

    # Show all the plots
    plt.show()

    
main()


