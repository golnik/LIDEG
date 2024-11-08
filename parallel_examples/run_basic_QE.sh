#!/bin/bash

#SBATCH -J graphene                # Job name
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks-per-node=1        # Number of tasks (CPUs) per node
#SBATCH --cpus-per-task=10          # Number of CPUs per task
#SBATCH --mem-per-cpu=5G           # Memory per CPU
#SBATCH --time=20:00:00             # Time limit (hh:mm:ss)
#SBATCH --account=graphene         # Account name
#SBATCH --partition=standard       # Partition name

######################################################

############  env  setting   please run this in puma
 
######################################################

# Load necessary modules
module purge
module load hdf5-intel/1.12.0 

# Add Quantum Espresso to the PATH
export PATH="/xdisk/ngolubev/mingruiyuan/quantum_espresso/q-e-qe-7.0/build/bin:$PATH"

# Save the current directory
ORIGINAL_PATH=$(pwd)

######################################################

############  Output path

######################################################

# Define the target path
TARGET_PATH="/xdisk/ngolubev/mingruiyuan/QE_diffr_time"

# Check if the target path exists, if not, create it
if [ ! -d "$TARGET_PATH" ]; then
    echo "Target path does not exist. Creating: $TARGET_PATH"
    mkdir -p "$TARGET_PATH"
fi

# Change to the target directory
cd "$TARGET_PATH"

# Clean up old directories and create necessary ones
mkdir ./wfc  ./outputs ./transition_re

# echo "Cleanup and setup completed at $TARGET_PATH."
# # Print environment variables for debugging
# echo "SLURM_NNODES = $SLURM_NNODES"
# echo "SLURM_NTASKS = $SLURM_NTASKS"
# echo "SLURM_NTASKS_PER_NODE = $SLURM_NTASKS_PER_NODE"

######################################################

############  QE_cal

######################################################

# Calculate the total number of processes (total CPUs used)
np=${SLURM_NTASKS:-$(($SLURM_NNODES * $SLURM_NTASKS_PER_NODE))}

# Check if np is set correctly
if [ -z "$np" ] || [ "$np" -le 0 ]; then
    echo "Error: process count should be greater than 0."
    exit 1
fi

# Run Quantum Espresso using mpirun
mpirun -np $np pw.x -inp "$ORIGINAL_PATH/inputs/graphene_scf.in" > outputs/graphene_scf.out

# Run additional Quantum Espresso tasks as needed
mpirun -np $np dos.x -inp "$ORIGINAL_PATH/inputs/graphene_dos.in" > outputs/graphene_dos.out
mpirun -np $np pw.x -inp "$ORIGINAL_PATH/inputs/graphene_bands.in" > outputs/graphene_bands.out
mpirun -np $np bands.x -inp "$ORIGINAL_PATH/inputs/graphene_bands_pp.in" > outputs/graphene_bands_pp.out

mpirun -np $np wfck2r.x < "$ORIGINAL_PATH/inputs/wfck2r.in" > outputs/wfck2r.out

###########################

module purge
module add intel/2020.4
module add gnu8/8.3.0
module add openmpi3/3.1.4
module add python/3.9/3.9.10

graphene_prog_path="/home/u18/mingruiyuan/LIDEG"

#input file
input=./inputs/dynamics_input.ini

$graphene_prog_path/build/Extract_QE_WF.exe $input  # get wavefunction from wfck2r.oct file psi_n(k,r)
$graphene_prog_path/build/Transition.exe $input     # calculate the Fourier transform of transition denisty for certain S

rm wfc* -rf

exit 0
