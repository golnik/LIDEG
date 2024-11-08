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
base_dir="/xdisk/ngolubev/mingruiyuan/QE_diffr_test"

# Nkx=65   #must be same with Nkx, Nky in inputs/dynamics_inout.ini and num_points in kpoints.py !!!!!!!!!!!!!!!!!!!!!
# Nkxy=$Nkx*$Nky

# Loop through t from 0 to 10
for t in {0..150}; do
    # Define the specific directory for this value of t
    target_dir="$base_dir/QE_diffr_$t"

    mkdir -p "$target_dir"
    
    # Copy ./inputs and ./field into the newly created folder
    cp -r ./inputs "$target_dir"
    cp -r ./field "$target_dir"
    # cp run_kspace.sh "$target_dir"
    cp run_basic_QE.sh "$target_dir"
    # cp plot_wfck2r.py "$target_dir"

    input_file_bands_pp="$target_dir/inputs/graphene_bands_pp.in"
    input_file_bands="$target_dir/inputs/graphene_bands.in"
    input_file_dos="$target_dir/inputs/graphene_dos.in"
    input_file_scf="$target_dir/inputs/graphene_scf.in"
    input_file_wfck2r="$target_dir/inputs/wfck2r.in"
    input_dynamics="$target_dir/inputs/dynamics_input.ini"
    file_basic_QE="$target_dir/run_basic_QE.sh"
    
    # Replace outdir and filband paths in graphene_bands_pp.in
    sed -i "s|outdir = .*|outdir = '$target_dir/tmp',|" "$input_file_bands_pp"
    sed -i "s|filband = .*|filband = '$target_dir/outputs/graphene_bands.dat',|" "$input_file_bands_pp"

    # Replace outdir and filband paths in graphene_band.in
    sed -i "s|outdir = .*|outdir = '$target_dir/tmp',|" "$input_file_bands"

    # Replace outdir and filband paths in graphene_dos.in
    sed -i "s|outdir = .*|outdir = '$target_dir/tmp',|" "$input_file_dos"
    sed -i "s|fildos = .*|fildos = '$target_dir/outputs/graphene_dos.dat',|" "$input_file_dos"

    # Replace outdir and filband paths in graphene_band.in
    sed -i "s|outdir = .*|outdir = '$target_dir/tmp',|" "$input_file_scf"

    # Replace outdir and filband paths in wfck2r.in
    sed -i "s|outdir = .*|outdir = '$target_dir/tmp',|" "$input_file_wfck2r"

    # Replace outdir and filband paths in wfck2r.in
    sed -i "s|TARGET_PATH=.*|TARGET_PATH=\"$target_dir\"|" "$file_basic_QE"

    # Replace outdir and filband paths in wfck2r.in
    sed -i "s|QE_path  = /xdisk/ngolubev/mingruiyuan/QE_diffr/|QE_path  = /xdisk/ngolubev/mingruiyuan/QE_diffr_time/QE_diffr_$t/|" "$input_dynamics"  

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

for t in {0..150}; do

    target_dir="$base_dir/QE_diffr_$t"
    cd "$target_dir"
    #puma
    sbatch run_basic_QE.sh                 # Run the QE script in the target directory

    cd "$ORIGINAL_PATH"

done

echo "Paths updated."
