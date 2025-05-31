import numpy as np
import matplotlib.pyplot as plt
import os

# Define the input file path
file_path = '/home/sriman/Documents/Thesis/Z_Matrix_elements.txt'
output_plot_path = '/home/sriman/Documents/Thesis/Z_vs_Mns_loglog.png'

# Check if the file exists
if not os.path.exists(file_path):
    print(f"Error: File not found at {file_path}")
    exit()

try:
    # Load data from the file, skipping the header row
    # Usecols specifies which columns to read (0-indexed)
    data = np.loadtxt(file_path, skiprows=1, usecols=(0, 1))

    # Separate Z and M_n_s values
    Z = data[:, 0]
    Mns = data[:, 1]

    # Ensure Mns values are positive for log plot by taking absolute value
    # Filter out any zero or negative values that might cause issues with log
    valid_indices = (Mns != 0) # Assuming Mns might be zero, though unlikely here
    Z_filtered = Z[valid_indices]
    Mns_abs_filtered = np.abs(Mns[valid_indices])

    if len(Z_filtered) == 0:
        print("Error: No valid (non-zero Mns) data points found for plotting.")
        exit()

    # Calculate the logarithms (using base 10 for interpretability)
    log_Z = np.log10(Z_filtered)
    log_Mns_abs = np.log10(Mns_abs_filtered)

    # Calculate the slope using linear regression (polyfit degree 1)
    # The slope is the first element of the returned array
    slope, intercept = np.polyfit(log_Z, log_Mns_abs, 1)

    print(f"Calculated slope of log10(Z) vs log10(|M_n_s|): {slope:.4f}")

    # Create the log-log plot
    plt.figure(figsize=(10, 6))
    plt.loglog(Z_filtered, Mns_abs_filtered, 'bo', label='Data points') # Use loglog directly

    # Plot the fitted line on the log-log scale
    # Generate Z values spanning the range for the line
    Z_line = np.logspace(np.log10(Z_filtered.min()), np.log10(Z_filtered.max()), 100)
    # Calculate corresponding Mns values using the log-log slope and intercept
    # log10(Mns) = slope * log10(Z) + intercept
    # Mns = 10**(slope * log10(Z) + intercept)
    Mns_line = 10**(slope * np.log10(Z_line) + intercept)
    plt.loglog(Z_line, Mns_line, 'r-', label=f'Fit (slope={slope:.2f})')

    # Add labels and title
    plt.xlabel('Atomic Number (Z)')
    plt.ylabel('Absolute Matrix Element |M_n_s|')
    plt.title('Log-Log Plot of Z vs |M_n_s|')
    plt.grid(True, which="both", ls="--")
    plt.legend()

    # Save the plot
    plt.savefig(output_plot_path)
    print(f"Plot saved to {output_plot_path}")

    # Optionally display the plot
    # plt.show()

except FileNotFoundError:
    print(f"Error: Could not find the file {file_path}")
except Exception as e:
    print(f"An error occurred: {e}")
