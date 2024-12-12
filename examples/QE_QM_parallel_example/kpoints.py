import sys

# Parse command-line arguments
t = int(sys.argv[1])
output_file = sys.argv[2]

# Define the number of points for each column
num_points = 129
Bx_0 = 1.56068822e+00 

# Generate linearly spaced values between 0 and 1 for the first and second columns
first_col_values = [i / (num_points - 1) for i in range(num_points)]
second_col_values = [i / (num_points - 1) for i in range(num_points)]

# Third column is always 0, fourth column is always 1
third_col = "0.0000000000"
fourth_col = "1"

# Read data from tfile.dat to get the replacement value
tfile = "output/tfile.dat"
with open(tfile, "r") as f:
    lines = f.readlines()

# Skip the header line and select the t-th row
if t + 1 < len(lines):
    selected_row = lines[t + 1].split()
    # Extract the value from the second column (Ax)
    ax_value = float(selected_row[1])
    ax_shfit = (ax_value/Bx_0)
else:
    raise ValueError(f"Row {t} is out of range in the tfile.dat.")

# Open the specified output file to write
with open(output_file, "w") as f:
    # Loop through all combinations of first and second columns
    for first_col in first_col_values:
        for second_col in second_col_values:
            # Format the output for each line, replacing 0.001*t with ax_value
            line = f"{first_col + ax_shfit:13.10f} {second_col + ax_shfit:13.10f} {third_col:13} {fourth_col}\n"
            f.write(line)

print(f"File {output_file} has been generated.")
