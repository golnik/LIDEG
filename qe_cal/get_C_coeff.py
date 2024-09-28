import h5py
import numpy as np

# Path to the HDF5 file containing wavefunctions
hdf5_file = "/tmp/mingrui_qe_cal_tmp/graphene.save/wfc1.hdf5"

# Open the output text file to save all the rearranged information
output_file = "./wf/wavefunction_with_G_vectors_info.txt"

with open(output_file, 'w') as file:
    # Open the HDF5 file in read mode
    with h5py.File(hdf5_file, 'r') as f:
        # Check if both MillerIndices and evc datasets are available
        if 'MillerIndices' in f.keys() and 'evc' in f.keys():
            miller_indices = f['MillerIndices'][:]
            evc = f['evc']

            # Number of bands (nbnd) and number of plane waves (npwx)
            nbnd, total_coeffs = evc.shape
            npwx = total_coeffs // 2  # Since each complex number is represented by 2 floats

            # Write total number of G vectors
            file.write(f"Total number of G vectors: {npwx}\n\n")

            # Loop over each band
            for band_idx in range(nbnd):
                file.write(f"band_index: {band_idx + 1}\n")
                file.write(f"{'G vector, [b1, b2, b3]':<30} {'C_G real':<30} {'C_G imag':<30}\n")

                # Loop over G-vectors for the current band
                for g in range(npwx):
                    # Extract the G vector from Miller indices
                    g_vector = miller_indices[g]

                    # Extract the real and imaginary parts of the wavefunction coefficient
                    real_part = evc[band_idx, 2 * g]
                    imag_part = evc[band_idx, 2 * g + 1]

                    # Write the G vector and the corresponding real/imag parts of the coefficient
                    file.write(f"{str(g_vector.tolist()):<30} {real_part:<30} {imag_part:<30}\n")

                file.write("\n")  # Separate the bands with a new line for readability
