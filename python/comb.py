import os
import sys
import argparse
import params

zones = [[[0, 0], [1, 0], [1, 1], [0, 1], [-1, 0], [-1, -1], [0, -1], 
          [2, 1], [1, 2], [-1, 1], [-2, -1], [-1, -2], [1, -1], 
          [2, 0], [2, 2], [0, 2], [-2, 0], [-2, -2], [0, -2]]]

def format_line(index, values):
    return ' '.join(f"{v:20.10f}" for v in [index] + values) + '\n'

def write_headers(file, zones):
    file.write(f"# {'Time[1]'.rjust(17)}")
    file.write(f"{'BZ0[2]'.rjust(18)}")
    file.write(f"{'BZ1[3]'.rjust(18)}")
    file.write(f"{'BZ2[4]'.rjust(18)}")
    file.write(f"{'BZ3[5]'.rjust(18)}")

    col_indx = 6
    for BZ in zones:
        for spot in BZ:
            m = str(spot[0])
            n = str(spot[1])
            label_str = f"BZ({m},{n})[{col_indx}]"
            file.write(label_str.rjust(18))
            col_indx += 1
    file.write("\n")

def parse_input_file(input_file):
    outdir = None
    Nt = None
    with open(input_file, 'r') as file:
        for line in file:
            if line.strip().startswith('outdir'):
                outdir = line.split('=')[1].strip()
            if line.strip().startswith('Nt'):
                Nt = int(line.split('=')[1].strip())
    return outdir, Nt

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', required=True)
    parser.add_argument('-output', required=False)

    args = parser.parse_args()
    ifname = args.input

    outdir, Nt = parse_input_file(ifname)

    if outdir is None or Nt is None:
        print("Required parameters (outdir, Nt) not found in the input file.")
        sys.exit(1)

    forder_name = args.output

    if os.path.exists(ifname):
        # Assuming params.analyze_input(ifname) initializes some parameters
        pass  # Replace this line with the actual function call if needed
    else:
        print("Input file does not exist!")
        sys.exit(1)

    # Set the base format for the filenames and the range
    filename_pattern = f"{outdir}/diff_{{:06d}}.dat"
    start_index = 1
    end_index = Nt

    # Name for the new file to store normalized data
    normalized_filename = f"{outdir}/normalized_data.dat"

    # Additional files for specific data
    filename_intra = f"{outdir}/diffr_intra.dat"
    filename_inter = f"{outdir}/diffr_inter.dat"
    filename_total = f"{outdir}/diffr_total.dat"

    # Open new files to write the divided data
    with open(filename_intra, 'w') as file_intra, \
         open(filename_inter, 'w') as file_inter, \
         open(filename_total, 'w') as file_total:

        # Write headers to each file
        write_headers(file_intra, zones)
        write_headers(file_inter, zones)
        write_headers(file_total, zones)

        # Iterate over the file numbers, read, normalize, and write the data in order
        for i in range(start_index, end_index + 1):
            filename = filename_pattern.format(i)
            try:
                with open(filename, 'r') as infile:
                    data_line = infile.readline().strip().split()
                    if not data_line:
                        # If the file is empty, write an empty line in the output files
                        file_intra.write(''.ljust(18))
                        file_inter.write('\n')
                        file_total.write('\n')
                    else:
                        data_values = [float(num) for num in data_line]
                        # Write divided data to each specific file
                        file_intra.write(format_line(data_values[0], [data_values[col] for col in range(1, len(data_values), 3)]))
                        file_inter.write(format_line(data_values[0], [data_values[col] for col in range(2, len(data_values), 3)]))
                        file_total.write(format_line(data_values[0], [data_values[col] for col in range(3, len(data_values), 3)]))

            except FileNotFoundError:
                print(f"File {filename} not found. Skipping.")
            except ValueError:
                print(f"Could not convert the data in {filename} to floats. Skipping.")

    print("Data normalization and division complete. Check respective files.")
