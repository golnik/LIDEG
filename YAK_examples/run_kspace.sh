#!/bin/bash

#SBATCH -J graphene
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00
#SBATCH --account=ltamp
#SBATCH --partition=standard

module purge
module add intel/2020.4
module add gnu8/8.3.0
module add openmpi3/3.1.4
module add python/3.9/3.9.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

graphene_prog_path="../"

#input file
input=input.ini

#extract number of steps
Nt=`cat $input | grep "Nt" | awk '{print $NF}'`

#reciprocal space calculations
$graphene_prog_path/build/main.exe $input 
$graphene_prog_path/build/analysis.exe $input

exit 0
