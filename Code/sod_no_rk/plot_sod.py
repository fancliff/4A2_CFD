#
#   plot_contours
#                               
#   Script to plot a converged flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *

# import os
# print(os.getcwd())

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
    g['v_ratio'] = g['v']/0.89095233

    # Specify the parameters to plot
    fieldnames = ['ro_ratio','v_ratio','p']; 
    colnames = ['Density Ratio','Velocity ratio','Pressure Ratio']

    j_mid = (g['nj'])//2
    # j_pos = 0
    j_pos = j_mid

    #Calculate and output the error metrics to the exact solution
    file_path = 'sod.raw'

    errors = calc_sod_error(g, file_path, j_pos)
    print()
    print(f'Density Error (ro): {errors["ro_error"]:.4e}')
    print(f'Pressure Error (p): {errors["p_error"]:.4e}\n')
    print(f'Total (1+2): {errors["ro_error"]+errors["p_error"]:.4e}')
    print(f'Total (1*2): {errors["ro_error"]*errors["p_error"]:.4e}\n')
    print(f'Added Pointwise Error: {errors["added_error"]:.4e}')
    print(f'Multiplied Pointwise Error: {errors["multiplied_error"]:.4e}\n')

    # Load the exact solution from the file
    exact_data = np.loadtxt(file_path, delimiter='\t', skiprows=1)
    x_exact, ro_exact, p_exact = exact_data[:, 0], exact_data[:, 1], exact_data[:, 2]
    vx_exact = exact_data[:,3]

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
            exact = ro_exact
        if name == 'v_ratio':
            t0 = np.zeros(g['ni'])
            exact = vx_exact/0.89095233
        if name == 'p':
            t0 = np.ones(g['ni'])
            t0[(g['ni']//2):g['ni']] *= 0.1
            exact = p_exact
        ax.plot(g['x'][:,j_pos],t0,label='t=0.0s',color='red')
        ax.plot(x_exact,exact,label='exact',color='green')

        # Add labels and title (fortran indexing for j_pos printout)
        ax.set_title(f'{colnames[n]} against x at j = {j_pos+1}', fontsize = 14)
        ax.set_xlabel('x', fontsize = 12)
        ax.set_ylabel(colnames[n], fontsize = 12)
        ax.legend(loc = 'best', fontsize = 12)

        ax.set_ylim(-0.05,1.25)
        ax.set_xlim(-0.05,1.05)

    # Show all the plots
    plt.show()
    
main()


