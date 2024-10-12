import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

# Define the path to the file
file_path = "../outputs/graphene_bands.dat"

# Initialize lists to store the K vectors and energy values
k_vectors = []
energy = []

# Open the file and read it line by line
with open(file_path, "r") as file:
    lines = file.readlines()

    # Skip the first line (header) and start processing from the second line
    for i in range(1, len(lines)):
        # Even-numbered lines (index 2, 4, 6, ...) are k-vectors
        if (i - 1) % 2 == 0:
            k_vector = lines[i].strip().split()
            k_vector = [float(value) for value in k_vector]
            k_vectors.append(k_vector)

        # Odd-numbered lines (index 1, 3, 5, ...) are energy values
        elif (i - 1) % 2 == 1:
            energy_values = lines[i].strip().split()
            energy_values = [float(value) for value in energy_values]
            energy.append(energy_values)

# Extract kx, ky from k_vectors and the 4th and 5th energy values from energy
kx = [k[0] for k in k_vectors]  # First component of k-vector (kx)
ky = [k[1] for k in k_vectors]  # Second component of k-vector (ky)
energy_4th = [e[3] for e in energy]  # 4th band energy values
energy_5th = [e[4] for e in energy]  # 5th band energy values

# Convert lists to numpy arrays
kx = np.array(kx)
ky = np.array(ky)
energy_4th = np.array(energy_4th)
energy_5th = np.array(energy_5th)

# Create a meshgrid for interpolation
kx_unique = np.linspace(min(kx), max(kx), 65)
ky_unique = np.linspace(min(ky), max(ky), 65)
kx_grid, ky_grid = np.meshgrid(kx_unique, ky_unique)

# Interpolate the energy values for the 4th and 5th bands onto the regular grid
energy_4th_grid = griddata((kx, ky), energy_4th, (kx_grid, ky_grid), method='cubic')
energy_5th_grid = griddata((kx, ky), energy_5th, (kx_grid, ky_grid), method='cubic')

# Create a 3D surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the 4th band surface
surf1 = ax.plot_surface(kx_grid, ky_grid, energy_4th_grid, cmap='viridis', alpha=0.7, label='4th band')

# Plot the 5th band surface with a different color map
surf2 = ax.plot_surface(kx_grid, ky_grid, energy_5th_grid, cmap='plasma', alpha=0.7, label='5th band')

# Set labels
ax.set_xlabel('kx')
ax.set_ylabel('ky')
ax.set_zlabel('Energy')

# Add a color bar for each surface
fig.colorbar(surf1, ax=ax, shrink=0.5, aspect=5, label='4th Band Energy')
fig.colorbar(surf2, ax=ax, shrink=0.5, aspect=5, label='5th Band Energy')

# Set the same aspect ratio for kx and ky by setting equal limits
ax.set_box_aspect([1, 1, 0.5])  # Aspect ratio: 1:1 for kx and ky, and shorter z-axis

# Manually set equal limits for kx and ky
kx_range = max(kx) - min(kx)
ky_range = max(ky) - min(ky)
max_range = max(kx_range, ky_range)

ax.set_xlim([min(kx), min(kx) + max_range])
ax.set_ylim([min(ky), min(ky) + max_range])

# Print min/max values for debugging purposes
print(f"Min of 5th band: {min(energy_5th)}")
print(f"Max of 4th band: {max(energy_4th)}")

# Show the plot
plt.show()
