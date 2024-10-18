import numpy as np
import os

# Define the file path to the .oct file
filepath = '../wfck2r.oct'

# Check if the file exists
if not os.path.isfile(filepath):
    raise FileNotFoundError(f'File not found: {filepath}')

# Open the .oct file
with open(filepath, 'r') as fid:

    # Initialize counters and flags
    isReadingUnkr = False
    wavefunctions = []
    dims = []

    # Loop through the file line by line
    for line in fid:
        line = line.strip()

        # Check if we are at the start of the unkr data block
        if '# name: unkr' in line:
            isReadingUnkr = True
            continue
        
        # Check if we've reached the dimensions block
        if isReadingUnkr and 'ndims:' in line:
            dims = list(map(int, fid.readline().split()))
            # dims will now contain [nr1x, nr2x, nr3x, nbands, nkpoints]
            continue
        
        # If we are reading the unkr matrix data, parse the complex values
        if isReadingUnkr and line.startswith('('):
            re, im = map(float, line[1:-1].split(','))  # Extract real and imaginary parts
            wavefunctions.append(complex(re, im))

# Convert wavefunctions list to a NumPy array
wavefunctions = np.array(wavefunctions)

# Extract dimensions: [nr1x, nr2x, nr3x, nbands, nkpoints]
nr1x, nr2x, nr3x, nbands, nkpoints = dims

# Ensure the number of wavefunctions matches the product of dimensions
expected_size = nr1x * nr2x * nr3x * nbands * nkpoints
if wavefunctions.size != expected_size:
    raise ValueError("Mismatch between wavefunction size and dimensions specified in the .oct file.")

# Reshape wavefunctions into a 5D matrix [nkpoints, nbands, nr3x, nr2x, nr1x]
wavefunctions = wavefunctions.reshape((nkpoints, nbands, nr3x, nr2x, nr1x))

# Ensure the output directory exists
output_dir = '../wf'
os.makedirs(output_dir, exist_ok=True)

# Loop over k-points and bands to save each wavefunction in a separate file
for k in range(nkpoints):
    for n in range(nbands):
        # Construct the filename for each wavefunction
        filename = os.path.join(output_dir, f'wfc_{n+1}_{k+1}.dat')
        
        # Extract the wavefunction for the current band and k-point
        wf = wavefunctions[k, n, :, :, :]
        
        # Flatten the wavefunction data
        wf_flat = wf.flatten()
        
        # Open the file for writing and save the real and imaginary parts of the wavefunction
        with open(filename, 'w') as fileID:
            for val in wf_flat:
                real_part = val.real
                imag_part = val.imag
                fileID.write(f'{real_part:.6e} {imag_part:.6e}\n')
        
        # Display a message for confirmation
        print(f'Saved wavefunction n={n+1} k={k+1} to {filename}')

print('All wavefunctions have been saved.')
