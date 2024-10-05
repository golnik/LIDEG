import os
import h5py
import numpy as np
import glob

# Define the directory containing HDF5 files and the output directory


hdf5_directory = "../tmp/graphene.save/"
output_directory = "../coeff"
os.makedirs(output_directory, exist_ok=True)

# Find all HDF5 files in the directory that match the pattern wfc*.hdf5
hdf5_files = glob.glob(os.path.join(hdf5_directory, "wfc*.hdf5"))

# Function to get all attributes of an HDF5 object
def get_attributes(hdf_object):
    attributes = {}
    for key, value in hdf_object.attrs.items():
        attributes[key] = value
    return attributes

# Loop over each HDF5 file found
for hdf5_file in hdf5_files:
    # Extract the base file name (without extension) to use for the output text file
    base_filename = os.path.splitext(os.path.basename(hdf5_file))[0]
    output_file = os.path.join(output_directory, f"coeff_{base_filename}.dat")

    with open(output_file, 'w') as file:
        # Open the HDF5 file in read mode
        with h5py.File(hdf5_file, 'r') as f:
            # Retrieve and write root attributes
            root_attributes = get_attributes(f)
            file.write("Root Attributes:\n")
            for key, value in root_attributes.items():
                file.write(f"{key}: {value}\n")
            file.write("\n")  # Add a newline for better readability

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

    #print(f"Processed {hdf5_file} and generated {output_file}")
