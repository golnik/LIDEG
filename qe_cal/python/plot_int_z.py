import os
import numpy as np
import matplotlib.pyplot as plt

# Initialize arrays to accumulate the sum of the real and imaginary parts
x_sum = None
y_sum = None
re_psi_integrated = None

# Loop over the files from 1 to 8
n = 5
m = 6

for band in range(n, m):
    file_path = f'../wf/rho_band_{band}.dat'
    
    # Load the data, skipping the header line
    data = np.loadtxt(file_path, skiprows=1)
    
    # Extract x, y, z, and the real part of the wavefunction (4th column)
    x = data[:, 0] 
    y = data[:, 1]
    z = data[:, 2]    
    re_psi = data[:, 3]
    
    # Perform integration over z (sum over z-axis)
    # Here, we assume that the data points are evenly spaced in z and
    # the z-axis is integrated by simple summation.
    # If z-spacing is uneven, use a more accurate integration method like np.trapz().
    
    # Create unique combinations of (x, y) grid points and sum over z
    xy_pairs, indices = np.unique(np.vstack([x, y]).T, axis=0, return_inverse=True)
    re_psi_2d = np.zeros_like(xy_pairs[:, 0])
    
    for idx in range(len(re_psi_2d)):
        re_psi_2d[idx] = re_psi[indices == idx].sum()  # Sum wavefunction over z for each (x, y)

    # Initialize the arrays on the first iteration
    if re_psi_integrated is None:
        x_sum = xy_pairs[:, 0]
        y_sum = xy_pairs[:, 1]
        re_psi_integrated = np.zeros_like(re_psi_2d)
    
    # Sum the integrated real part of the wavefunction over the bands
    re_psi_integrated += re_psi_2d

# Create a scatter plot with color mapping for the summed real part
plt.figure(figsize=(8, 6))
plt.tricontourf(x_sum, y_sum, re_psi_integrated, levels=100, cmap='coolwarm')  # Smooth the data with contour
plt.colorbar(label='Summed Wavefunction Amplitude (Integrated over z)')
plt.title(f'Summed 2D Wavefunction for Graphene (Bands {n} to {m}, Integrated over z)')
plt.xlabel('X')
plt.ylabel('Y')
plt.gca().set_aspect('equal', adjustable='box')

# Save the plot
os.makedirs('../figures', exist_ok=True)
plt.savefig(f'../figures/{n}_{m}_intz_wavefunction_integrated.png')

# Show the plot (uncomment if you want to display the plot interactively)
# plt.show()
