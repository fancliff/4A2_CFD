#
#   plot_spacetime
#                               
#   Script to plot a space-time-mach plot on the centerline from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_contours.py casename"

# Import modules and functions
from routines import *
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from routines import read_case, calc_secondary

def load_mach_frames(file_list, av):
    """
    Load only Mach number data for all frames from the provided list of files.
    Returns a list of Mach arrays
    """
    mach_data_list = []
    nj = None

    for file in file_list:
        # Read the data for the current frame
        g = read_case(file)
        g = calc_secondary(av, g)  # Calculate secondary variables like Mach

        # Store only the Mach field
        mach_data_list.append(g['mach'])

        if nj is None:
            # Calculate centerline index
            nj = g['nj']
            centerline_idx = (nj + 1) // 2 - 1  # Zero-based centerline index
            x_coords = g['x'][:,centerline_idx]

    return mach_data_list, x_coords, centerline_idx

def generate_space_time_mach_plot():
    """
    Generates a space-time Mach plot for the centerline using preloaded Mach data.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import glob

    # File handling
    folder = 'tunnel_blowdown'
    file_pattern = os.path.join(folder, 'out_unste_*.bin')
    files = sorted(glob.glob(file_pattern))

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    av = read_settings(inname)
    
    if not files:
        print(f"No files found matching {file_pattern}.")
        return
    
    print(f"Found {len(files)} files. Preloading Mach data...")

    # User inputs
    num_total_steps = int(input('Enter Total Timesteps: '))  # Total number of steps in the solution
    timestep = 1.63889786E-06       # Time step per solution step
    num_frames = len(files)         # Number of frames to process
    
    # Load Mach data, x_coords, and centerline index
    mach_data, x_coords, centerline_idx = load_mach_frames(files, av)
    
    # Time array for the frames
    frame_times = np.linspace(0, num_total_steps * timestep, num_frames)

    # Print shape of frame_times and x_coords
    print(f"Shape of frame_times: {frame_times.shape}")
    print(f"Shape of x_coords: {x_coords.shape}")

    # Extract Mach data along the centerline
    centerline_mach = []
    for mach in mach_data:
        centerline_mach.append(mach[:, centerline_idx])  # Fixed j, varying i
    
    # Convert to a 2D NumPy array
    centerline_mach = np.array(centerline_mach)

    # Print shape of centerline Mach data
    print(f"Shape of centerline_mach: {centerline_mach.shape}")

    # Generate meshgrid for x (spatial) and time
    X, T = np.meshgrid(x_coords, frame_times)

    # Print shape of the meshgrid
    print(f"Shape of X: {X.shape}")
    print(f"Shape of T: {T.shape}")

    # Check alignment of dimensions
    if X.shape != centerline_mach.shape or T.shape != centerline_mach.shape:
        print("Error: Meshgrid dimensions do not match the Mach data dimensions.")
        return
    else:
        print("Meshgrid and Mach data dimensions match. Proceeding with plot generation.")

    # Generate the space-time Mach plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Use pcolormesh with the corrected grid
    c = ax.pcolormesh(X, T, centerline_mach, shading='auto', cmap='viridis', vmin=0, vmax=3)

    # Add labels and colorbar
    fig.colorbar(c, ax=ax, label='Mach Number')
    ax.set_xlabel('Spatial Coordinate (x)')
    ax.set_ylabel('Time (s)')
    ax.set_title('Space-Time Mach Plot Along the Centerline')

    # Save and show the plot
    plt.tight_layout()
    plt.savefig('space_time_mach_plot.png', dpi=150)
    plt.show()


# Run the function
generate_space_time_mach_plot()
