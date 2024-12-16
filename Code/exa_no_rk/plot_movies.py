#
#   plot_movies
#                               
#   Script to plot a unsteady movie flowfield from the 4A2 solver
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/plot_movies.py casename"

# Import modules and functions
import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from routines import *

# Function to process files and generate movies
def generate_movies():

    # Construct full filenames to read the run data
    inname = 'input_' + sys.argv[-1] + '.txt'
    av = read_settings(inname)

    # Reference pressures are calculated only once and use the final values
    outname = 'out_final_' + sys.argv[-1] + '.bin'
    g = read_case(outname)
    g = calc_secondary(av, g)
    g_inlet = cut_i(g, 0)
    pstag_ref, _ = mass_av(g_inlet, 'pstag')
    p_ref, _ = area_av(g_inlet, 'p')
    print('')
    print(f'Reference static pressure: {p_ref:.0f}, Reference stagnation pressure: {pstag_ref:.0f}\n')

    # Path to the directory containing the input files
    folder = 'tunnel_blowdown'
    file_pattern = os.path.join(folder, 'out_unste_*.bin')
    files = sorted(glob.glob(file_pattern))

    if not files:
        print(f"No files found in the directory {folder} matching the pattern {file_pattern}.")
        return

    print(f"Found {len(files)} files to process.")

    # Parameters to plot
    fieldnames = ['mach', 'cp', 'cpstag']
    colnames = ['Mach number', 'Static pressure coefficient', 'Stagnation pressure coefficient']

    # Only preload every n-th frame
    n = 1  # Modify this to control how many frames to skip

    # Preload every n-th frame
    for field_idx, (current_field, current_colname) in enumerate(zip(fieldnames, colnames)):

        frames_data = []
        print(f"Preloading frames for {current_field}...")

        for idx, file in enumerate(files):
            if idx % n == 0:  # Only load every n-th frame
                g = read_case(file)

                # Calculate secondary variables
                g = calc_secondary(av, g)

                g['cp'] = (g['p'] - p_ref) / (pstag_ref)
                g['cpstag'] = (g['pstag'] - pstag_ref) / (pstag_ref)

                frames_data.append(g)

        print(f"Preloaded {len(frames_data)} frames for {current_field}.")

        # Create a figure for animation
        fig, ax = plt.subplots(figsize=[9.6, 7.2], dpi=80)
        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        # Plot the walls only once (static)
        plot_wall(ax, g)

        # Placeholder for the color mesh and contour lines
        hc = None
        contour_lines = None

        # Function to initialize the animation
        def init():
            ax.set_aspect('equal', adjustable='box')
            ax.axis('off')
            return []

        # Function to update each frame
        def update(frame_idx):
            nonlocal hc, contour_lines

            # Print a counter so that progress can be monitored
            print(frame_idx)

            # Get the preloaded frame
            g = frames_data[frame_idx]

            # Determine vmin and vmax based on the field
            if current_field == 'cp':
                vmin, vmax = -2, 2
            elif current_field == 'mach':
                vmin, vmax = 0, 3  # Typical Mach number range
            elif current_field == 'cpstag':
                vmin, vmax = -2, 2
            else:
                vmin, vmax = None, None  # Auto-scale if undefined

            # Clear the axis for dynamic content, but not the walls
            ax.clear()
            ax.set_aspect('equal', adjustable='box')
            ax.axis('off')

            # Plot the color mesh with fixed scale
            hc = ax.pcolormesh(
                g['x'], g['y'], g[current_field], shading='gouraud', vmin=vmin, vmax=vmax
            )
            colorbar(hc, current_colname)
            
            # Add Mach = 1 contour only if levels exist in the range
            if current_field == 'mach' and np.min(g['mach']) <= 1.0 <= np.max(g['mach']):
                ax.contour(g['x'], g['y'], g['mach'], [1.0], colors='w', linewidths=0.5)

                
            '''
            # Add cpstag contours only if data range supports it
            if current_field == 'cpstag' and np.min(g['cpstag']) < 2 and np.max(g['cpstag']) > -2:
                levels = np.arange(-2, 2, 0.02)
                contour_lines = ax.contour(g['x'], g['y'], g['cpstag'], levels=levels, colors='w', linewidths=0.5)
                ax.clabel(contour_lines, inline=False, fontsize=8)
            '''

            return [hc, contour_lines] if contour_lines else [hc]

        # Number of frames is the total number of preloaded frames
        num_frames = len(frames_data)

        ani = FuncAnimation(fig, update, frames=num_frames, init_func=init, blit=False)

        try:
            output_video = f'tunnel_blowdown_movie_{current_field}.mp4'

            # Adjust figure size and layout
            fig.set_size_inches(8, 6)  # Set to reasonable dimensions
            fig.tight_layout()  # Prevent label overlap

            ani.save(output_video, fps=10, extra_args=['-vcodec', 'mpeg4'])
            print(f"Movie for {current_field} saved as {output_video}")
        except Exception as e:
            print(f"Failed to save {current_field} as MP4 due to: {e}")
            print("Falling back to GIF format...")

            # Reduce resolution for GIF saving
            fig.set_size_inches(6.4, 4.8)  # Smaller figure size
            ani.save(f'tunnel_blowdown_movie_{current_field}.gif', fps=10, writer='pillow')
            print(f"Movie for {current_field} saved as tunnel_blowdown_movie_{current_field}.gif")

        # Close the figure to free memory
        plt.close(fig)

# Call the function
generate_movies()






