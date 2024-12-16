#
#   generate_bcs                      
#                               
#   Script to generate time varying boundary conditions
#
#   Change to the directory you want to execute the script within and execute 
#   with "python path_to_script/generate_bcs.py"

# Import modules and functions
import numpy as np
from routines import *

import numpy as np

def generate_bcs(initial_pressure,
                 total_time,
                 decay_time,
                 num_points,
                 output_file,
                 atmospheric_pressure=100000,
                 delay_time=0.1,
                 hold_time=5.0
):
    """
    Generates a .txt file with time, stagnation pressure values (constant during hold time, then decaying),
    and a constant atmospheric pressure column.

    Parameters:
        initial_pressure (float): Initial stagnation pressure (P0_initial).
        total_time (float): Total time sampled (in seconds).
        decay_time (float): Time for stagnation gauge pressure to reach 1%
        num_points (int): Number of data points to generate.
        output_file (str): Name of the output .txt file.
        atmospheric_pressure (float): Constant atmospheric pressure to add (default: 100,000 Pa).
        delay_time (float): Time delay before the valve opens (default: 0.1 s).
        hold_time (float): Duration during which stagnation pressure remains constant before decay (default: 5 s).
    """
    # Calculate the decay constant K to ensure P0_final = 0.01 * P0_initial after total_time
    K = -np.log(0.01) / decay_time  # Decay constant for 1% decay after decay_time

    # Generate time points from 0 to total_time
    time = np.linspace(0, total_time, num_points)

    # Initialize stagnation pressure array
    stagnation_pressure = np.zeros(num_points)

    # Set stagnation pressure to atmospheric pressure for the first timestep
    stagnation_pressure[0] = atmospheric_pressure

    # Compute stagnation pressure for times after the valve is opened
    for i in range(1, num_points):
        if time[i] < delay_time:
            stagnation_pressure[i] = atmospheric_pressure
        elif time[i] < delay_time + hold_time:
            # During the holding period, the stagnation pressure remains constant at the initial pressure
            # This can be thought of as the effect of a pressure regulator
            stagnation_pressure[i] = atmospheric_pressure + initial_pressure
        else:
            # After hold_time, the stagnation pressure starts decaying
            stagnation_pressure[i] = initial_pressure * np.exp(-K * (time[i] - delay_time - hold_time)) + atmospheric_pressure

    # Write to the output file
    with open(output_file, "w") as file:
        for t, p in zip(time, stagnation_pressure):
            # Write time, stagnation pressure, and constant atmospheric pressure
            file.write(f"{t:.6f}    {p:.6f}    {atmospheric_pressure:.6f}\n")

    print(f"File '{output_file}' successfully created with {num_points} points.")


# Main script logic
if __name__ == "__main__":
    # User-defined parameters
    initial_pressure = float(input("Enter initial stagnation GAUGE pressure (P0_initial): "))
    total_time = 40.0  # Total time sampled in seconds
    decay_time = 30.0  # Time for decay to 1% of gauge pstag
    num_points = int(input("Enter number of points to generate: "))
    output_file = 'bcs_tunnel.txt'
    # output_file = input("Enter output file name (e.g., 'bcs.txt'): ")
    atmospheric_pressure = 100000
    delay_time = float(input("Enter valve opening delay time (in seconds): "))
    hold_time = float(input("Enter holding time (in seconds): "))

    # Generate the file
    generate_bcs(initial_pressure,
                 total_time,
                 decay_time,
                 num_points,
                 output_file,
                 atmospheric_pressure,
                 delay_time,
                 hold_time)



