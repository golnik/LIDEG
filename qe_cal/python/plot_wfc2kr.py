import numpy as np
import matplotlib.pyplot as plt
import os

# Unit cell parameters
a = 2.46
A1 = np.array([a*2/np.sqrt(3),  a/np.sqrt(3),  0.0])  # first lattice vector
A2 = np.array([a*2/np.sqrt(3), -a/np.sqrt(3),  0.0]) # second lattice vector
A3 = np.array([0.0,   0.0,  10.0])  # third lattice vector (for z-direction)

# Grid parameters
nr1x = 15  # grid size in x direction
nr2x = 15  # grid size in y direction
nr3x = 54  # grid size in z direction
nkpoints = 33 # number of k-points
nbands = 8  # number of bands
n = 4  # select the band for which you want the 2D plot

# Initialize rho_n(r) to zero
rho   = np.zeros((nr1x, nr2x))
rho_n = np.zeros((nr1x, nr2x))

# for sum_n in range(1, n+1):
for k in range(1, nkpoints + 1):
    # Load the wavefunction data for this k-point and band n
    file_name = f'../wf/wfc_{n}_{k}.dat'
    if not os.path.isfile(file_name):
        raise FileNotFoundError(f"File not found: {file_name}")

    data = np.loadtxt(file_name)  # assuming the file has two columns: real and imaginary parts

    # Reconstruct the complex wavefunction and reshape to (nr3x, nr2x, nr1x)
    psi_n_k = data[:, 0] + 1j * data[:, 1]
    psi_n_k = psi_n_k.reshape((nr3x, nr2x, nr1x))

    norm = np.sqrt(np.sum(np.abs(psi_n_k)**2))
    
    # Normalize the wavefunction
    # psi_n_k /= norm
    
    # Print the value for verification (optional)
    # print(psi_n_k[0, 0, 0])

    # Integrate over the z-direction and sum the charge density
    for x in range(nr1x):
        for y in range(nr2x):
            rho_n[x, y] += np.sum(np.abs(psi_n_k[:, y, x])**2)
        
        #rho_n /= nkpoints

    # print(np.sqrt(np.sum(rho_n)))

    # rho = rho + rho_n


# Normalize the charge density (optional)


# Calculate real space coordinates
X, Y = np.meshgrid(np.arange(nr1x), np.arange(nr2x))  # indices for grid
realX = X * A1[0] / nr1x + Y * A2[0] / nr2x
realY = X * A1[1] / nr1x + Y * A2[1] / nr2x

# Plot the 2D contour of the charge density in real space
plt.figure()
contour = plt.contourf(realX, realY, rho_n, levels=20)  # 20 levels of contours
plt.colorbar(contour)
plt.xlabel('x (Å)')
plt.ylabel('y (Å)')
plt.title(f'2D Contour Plot of Charge Density for Band {n}')
plt.show()
