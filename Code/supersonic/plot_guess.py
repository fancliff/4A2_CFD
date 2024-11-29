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
    inname = 'input_' + sys.argv[-1] + '.txt'
    filename = 'out_guess_' + sys.argv[-1] + '.bin'

    # Read the case from file
    g = read_case(filename)
    av = read_settings(inname)
    g = calc_secondary(av,g)

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
 
        # Plot filled contour levels
        hc = ax[n].pcolormesh(g['x'],g['y'],g[name],shading='gouraud')

  	# Add colorbar with variable name
        colorbar(hc,name)

        # Draw the walls of the block
        plot_wall(ax[n],g)


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

        # Add Mach = 1 contours
        if name == 'mach':
            ax.contour(g['x'],g['y'],g['mach'],[1.0],colors='w',
                linewidths=0.5)
            
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


