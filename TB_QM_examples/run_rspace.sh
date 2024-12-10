#!/bin/bash

#SBATCH -J graphene
#SBATCH --nodes=22
#SBATCH --ntasks-per-node=7
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=152
#SBATCH --mem-per-cpu=4G
#SBATCH --time=40:00:00
#SBATCH --account=ltamp
#SBATCH --partition=standard

module purge
module add gnu8/8.3.0
module add openmpi3/3.1.4
module add ohpc
module add python/3.9/3.9.10

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

graphene_prog_path="../"

#input file
input=input.ini

#extract number of steps
Nt=`cat $input | grep "Nt" | awk '{print $NF}'`

#create MPI tasks list
MPI_tasks_list_fname="$PWD/MPI_tasks_list"
echo "$graphene_prog_path/build/data_writer.exe -f $PWD/$input --rspace" > $MPI_tasks_list_fname
for (( i=1; i<=$Nt; i++ )); do
    echo "$graphene_prog_path/build/rspace.exe $PWD/$input $i" >> $MPI_tasks_list_fname
done

MPI_script="$graphene_prog_path/mpi_manager/manager.py"
outdir="$PWD/MPI_out"

#start real space calculations
mpiexec python3 $MPI_script $MPI_tasks_list_fname $outdir

exit 0
