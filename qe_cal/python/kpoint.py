# Define the number of points for each column
num_points = 13

# Generate linearly spaced values between 0 and 1 for the first and second columns
first_col_values = [i / (num_points - 1) for i in range(num_points)]
second_col_values = [i / (num_points - 1) for i in range(num_points)]

# Third column is always 0, fourth column is always 1
third_col = "0.0000000000"
fourth_col = "1"

# Open a file to write the output
with open("output.txt", "w") as f:
    # Loop through all combinations of first and second columns
    for second_col in second_col_values:
        for first_col in first_col_values:
            # Format the output for each line
            line = f"{first_col:13.10f} {second_col:13.10f} {third_col:13} {fourth_col}\n"
            f.write(line)

print("File output.txt has been generated.")
