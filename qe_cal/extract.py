import os
import h5py
import numpy as np
import re


# Step 1: Extract k-point and weight information from the SCF output file
def extract_kpoints_from_scf_output(file_path):
    kpoints = []
    pattern = r"k\(\s*\d+\)\s*=\s*\(\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*\),\s*wk\s*=\s*(\d+\.\d+)"
    
    with open(file_path, 'r') as f:
        for line in f:
            match = re.search(pattern, line)
            if match:
                # Extract kx, ky, kz, and weight as floats
                kx, ky, kz, weight = map(float, match.groups())
                # Append the k-point to the list
                kpoints.append({'k': (kx, ky, kz), 'weight': weight})

    # Determine the midpoint of the kpoints list
    half_index = len(kpoints) // 2

    # Slice the list to keep only the second half
    kpoints = kpoints[half_index:]

    # Save the filtered kpoints to a file in the "wf" folder
    output_path = 'wf/kpoints'
    with open(output_path, 'w') as f:
        for point in kpoints:
            kx, ky, kz = point['k']
            weight = point['weight']
            f.write(f"{kx:.6f} {ky:.6f} {kz:.6f} {weight:.6f}\n")

    return kpoints

# Step 2: Compute the wavefunction for each k-point at a given real-space point
def compute_wavefunction_at_r(k_vec, g_vectors, coeff_real, coeff_imag, r_vec):
    wavefunction = 0.0 + 0.0j  # Initialize as a complex number

    # Loop through all G-vectors and their corresponding coefficients
    tmp = 0
    for g, real, imag in zip(g_vectors, coeff_real, coeff_imag):
        tmp = tmp + 1
        
        # Assuming g is a 3D vector (m1, m2, m3)
        m1, m2, m3 = g  # Decompose the g vector into components

        #print(g)
        
        # Apply the transformation: (m1+m2)/sqrt(3), (m1-m2), 0
        g_new = np.array([(m1 + m2) / np.sqrt(3), m1 - m2, 0])

        #print(g_new)
        
        g_plus_k = k_vec + g_new  # Calculate G + k with the transformed G vector
        phase_factor = np.exp(1j * np.dot(g_plus_k, (r_vec)) * 2 * np.pi)  # e^{i(k+G)·r}
        
        coefficient = real + 1j * imag  # Complex coefficient C_G
        wavefunction += coefficient * phase_factor  # Sum over the plane waves
        # print(tmp)
        # break
    
    return wavefunction

# Step 3: Summing over k-points to get the total wavefunction psi_n(r)
def sum_wavefunctions_over_kpoints(kpoints_data, hdf5_directory, real_space_grid, nbnd):
    total_wavefunction = np.zeros((len(real_space_grid), nbnd), dtype=complex)
    
    tmp = 0
    for i, kp in enumerate(kpoints_data):
        tmp = tmp + 1
        # print(tmp)
        # Construct the filename for the current k-point
        wfc_file = f"{hdf5_directory}/wfc{i+1}.hdf5"

        # Open the wavefunction HDF5 file
        with h5py.File(wfc_file, 'r') as f:
            wavefunctions = f['evc'][:]  # Assuming 'evc' contains wavefunctions
            
            # Extract the G-vectors from MillerIndices
            g_vectors = f['MillerIndices'][:]  # Reciprocal lattice vectors (G-vectors)

            # Loop over each band and sum wavefunctions weighted by the k-point weight
            for band_idx in range(nbnd):
                for r_idx, r_vec in enumerate(real_space_grid):
                    # Compute the wavefunction for this k-point
                    psi_k = compute_wavefunction_at_r(
                        np.array(kp['k']), g_vectors, 
                        wavefunctions[band_idx, 0::2], wavefunctions[band_idx, 1::2], 
                        r_vec
                    )
                    # Accumulate the weighted wavefunction
                    total_wavefunction[r_idx, band_idx] += kp['weight']/2 * psi_k
                    # kp['weight']/2 * psi_k

    return total_wavefunction

# Step 4: Save the wavefunction for each band and real-space point to the folder 'wf'
def save_wavefunction_to_files(total_wavefunction, real_space_grid, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Loop over each band and save wavefunction to a file
    for band_idx in range(total_wavefunction.shape[1]):
        filename = os.path.join(output_folder, f"wavefunction_band_{band_idx+1}.txt")
        with open(filename, 'w') as wf_file:
            wf_file.write(f"{'x':<12} {'y':<12}  {'Real part':<20} {'Imaginary part':<20} {'Norm part':<20}\n" ) 
            
            for r_idx, r_vec in enumerate(real_space_grid):
                wf_file.write(f"{r_vec[0]:<12.6f} {r_vec[1]:<12.6f} "
                              f"{total_wavefunction[r_idx, band_idx].real:<20.10f} "
                              f"{total_wavefunction[r_idx, band_idx].imag:<20.10f}"
                              f"{np.sqrt((total_wavefunction[r_idx, band_idx].real)**2 + (total_wavefunction[r_idx, band_idx].imag)**2):<20.10f}\n")

    print(f"Wavefunctions saved in {output_folder}")

# Step 5: Main procedure to read k-points, read wavefunctions and sum over k-points
def main(scf_output_file, hdf5_directory, real_space_grid, nbnd, output_folder):
    # Step 1: Extract k-point information and weights
    kpoints_data = extract_kpoints_from_scf_output(scf_output_file)

    #print(kpoints_data)
    # Step 2 & 3: Sum the wavefunctions over all k-points
    total_wavefunction = sum_wavefunctions_over_kpoints(kpoints_data, hdf5_directory, real_space_grid, nbnd)

    # Step 4: Save the total wavefunction to files
    save_wavefunction_to_files(total_wavefunction, real_space_grid, output_folder)

# Example usage
scf_output_file = './outputs/graphene_scf.out'  # Path to the SCF output file
hdf5_directory = '/tmp/mingrui_qe_cal_tmp/graphene.save/'  # Path to directory containing wfcXXX.hdf5 files

# Define the real-space grid (e.g., based on a 16x16 grid)
N = 4

Nx, Ny = 16*N, 16*N
real_space_grid = []
for i in range(Nx):
    for j in range(Ny):
        frac_coords = np.array([i / Nx, j / Ny, 0.0])  # Example grid in fractional coordinates
        # Assuming a unit cell to convert fractional coordinates to real-space
        cell_parameters = np.array([[N/np.sqrt(3), N/2, 0.0], [N/np.sqrt(3), -N/2, 0.0], [0.0, 0.0, 0.0]])  # unit of a = 2.46
        # cell_parameters = np.array([[2.13, 1.23, 0.0], [2.13, -1.23, 0.0], [0.0, 0.0, 10.0]])  # x y z
        real_space_grid.append(np.dot(frac_coords, cell_parameters))

# Set the number of bands
nbnd = 8  # Adjust based on your system

# Folder to save wavefunction files
output_folder = "wf"

# Run the main function
main(scf_output_file, hdf5_directory, real_space_grid, nbnd, output_folder)
