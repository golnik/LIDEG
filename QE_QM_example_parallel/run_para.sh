#!/bin/bash

#SBATCH -J graphene                # Job name
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks-per-node=1        # Number of tasks (CPUs) per node
#SBATCH --cpus-per-task=2          # Number of CPUs per task
#SBATCH --mem-per-cpu=5G           # Memory per CPU
#SBATCH --time=1:00:00             # Time limit (hh:mm:ss)
#SBATCH --account=graphene         # Account name
#SBATCH --partition=standard       # Partition name


# Set the base directory

ORIGINAL_PATH=$(pwd)
base_dir="/xdisk/ngolubev/mingruiyuan/QE_diffr_time"

# Loop through t from 0 to 10
for t in {0..301}; do
    # Define the specific directory for this value of t
    target_dir="$base_dir/QE_diffr_$t"

    mkdir -p "$target_dir"
    
    # Copy ./inputs and ./field into the newly created folder
    cp -r ./inputs "$target_dir"
    # cp -r ./field "$target_dir"
    # cp run_kspace.sh "$target_dir"
    cp run_basic_QE.sh "$target_dir"
    # cp plot_wfck2r.py "$target_dir"

    # input_file_bands_pp="$target_dir/inputs/graphene_bands_pp.in"
    input_file_bands="$target_dir/inputs/graphene_bands.in"
    # input_file_dos="$target_dir/inputs/graphene_dos.in"
    input_file_scf="$target_dir/inputs/graphene_scf.in"

    # Run the Python script to generate output.txt in the current QE_diffr_t directory
    python3 kpoints.py "$t" "$target_dir/output.txt"

    # Append the content of output.txt to graphene_scf.in
    cat "$target_dir/output.txt" >> "$input_file_bands"

    # Append the content of output.txt to graphene_scf.in
    cat "$target_dir/output.txt" >> "$input_file_scf"

    # cd "$target_dir"
    # #puma
    # sbatch run_basic_QE.sh                 # Run the QE script in the target directory

    # cd "$ORIGINAL_PATH"
done

for t in {0..301}; do
    target_dir="$base_dir/QE_diffr_$t"
    cd "$target_dir"
    #puma
    sbatch run_basic_QE.sh                 # Run the QE script in the target directory

    cd "$ORIGINAL_PATH"

done

echo "Paths updated."
