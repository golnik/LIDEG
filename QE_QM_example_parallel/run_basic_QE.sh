#!/bin/bash

#SBATCH -J graphene                # Job name
#SBATCH --nodes=1                  # Number of nodes
#SBATCH --ntasks-per-node=1        # Number of tasks (CPUs) per node
#SBATCH --cpus-per-task=50          # Number of CPUs per task
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

############  Calculation/tmp path

######################################################

# Define the target tmp path

TMP_DIR="/tmp/$SLURM_JOB_ID"
mkdir -p $TMP_DIR
echo "$TMP_DIR"

cp -r ./inputs "$TMP_DIR"
# cp -r ./field "$target_dir"
# cp run_kspace.sh "$target_dir"
# cp run_basic_QE.sh "$TMP_DIR"

# Change to the target directory
cd "$TMP_DIR"

# Create necessary folders
mkdir ./wfc  ./outputs ./transition_re

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
mpirun -np $np pw.x -inp "./inputs/graphene_scf.in" > outputs/graphene_scf.out

mpirun -np $np dos.x -inp "./inputs/graphene_dos.in" > outputs/graphene_dos.out

mpirun -np $np pw.x -inp "./inputs/graphene_bands.in" > outputs/graphene_bands.out

mpirun -np $np bands.x -inp "./inputs/graphene_bands_pp.in" > outputs/graphene_bands_pp.out

mpirun -np $np wfck2r.x < "./inputs/wfck2r.in" > outputs/wfck2r.out

######################################################

############  Extract wavefunction from wfck2r.oct and calculate transition

######################################################

module purge
module add intel/2020.4
module add gnu8/8.3.0
module add openmpi3/3.1.4
module add python/3.9/3.9.10

# Extract&Transition code path
graphene_prog_path="/home/u18/mingruiyuan/LIDEG_ALL/QE_QM_Diffr"

#input file
input=./inputs/dynamics_input.ini

$graphene_prog_path/build/Extract_QE_WF.exe $input  # get wavefunction from wfck2r.oct file psi_n(k,r)
$graphene_prog_path/build/Transition.exe $input     # calculate the Fourier transform of transition denisty for certain S

# copy the transition back to origional path
cp -r ./transition_re "$ORIGINAL_PATH"

exit 0
