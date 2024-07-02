
import os
# Set the base format for the filenames and the range
filename_pattern = "output/diff_{:06d}.dat"
start_index = 1
end_index = 151

# Name for the new file to store normalized data
normalized_filename = "output/normalized_data.dat"

# Additional files for specific data
filename_intra = "output/diffr_intra.dat"
filename_inter = "output/diffr_inter.dat"
filename_total = "output/diffr_total.dat"

# Read the first file to obtain normalization factors for columns a[2] to a[6]
try:
    with open(filename_pattern.format(start_index), 'r') as first_file:
        # Read the first line and split it into a list of floats
        normalization_factors = [float(num) for num in first_file.readline().strip().split()]
except FileNotFoundError:
    print(f"File {filename_pattern.format(start_index)} not found. Please check that the file exists.")
    raise
except ValueError:
    print(f"Could not convert the line of {filename_pattern.format(start_index)} to floats.")
    raise

# Define a function for consistent spacing in file writing
def format_line(index, values):
    return ' '.join(f"{v:20.6f}" for v in [index] + values) + '\n'

# Open new files to write the divided data
with open(normalized_filename, 'w') as outfile, \
     open(filename_intra, 'w') as file_intra, \
     open(filename_inter, 'w') as file_inter, \
     open(filename_total, 'w') as file_total:

    # Iterate over the file numbers, read, normalize, and write the data in order
    for i in range(start_index, end_index + 1):
        filename = filename_pattern.format(i)
        try:
            with open(filename, 'r') as infile:
                data_line = infile.readline().strip().split()
                data_values = [float(num) for num in data_line]
                # Normalize data based on the normalization factors
                #normalized_values = [data / factor for data, factor in zip(data_values, normalization_factors)]
                outfile.write(format_line(data_values[0], data_values))

                # Write divided data to each specific file
                file_intra.write(format_line(data_values[0], [data_values[col] for col in range(1, len(data_values), 3)]))
                file_inter.write(format_line(data_values[0], [data_values[col] for col in range(2, len(data_values), 3)]))
                file_total.write(format_line(data_values[0], [data_values[col] for col in range(3, len(data_values), 3)]))

        except FileNotFoundError:
            print(f"File {filename} not found. Skipping.")
        except ValueError:
            print(f"Could not convert the data in {filename} to floats. Skipping.")

# Remove all original diff_{:06d}.dat files
#for i in range(start_index, end_index + 1):
#    filename = filename_pattern.format(i)
#    try:
#        os.remove(filename)
#        #print(f"Removed file: {filename}")
#    except FileNotFoundError:
#        print(f"File {filename} not found when trying to delete. It may have already been deleted.")


print("Data normalization and division complete. Check respective files.")
