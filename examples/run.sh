#!/bin/bash

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

graphene_prog_path="../"

#input file
input=input.ini

#extract number of steps
Nt=`cat $input | grep "Nt" | awk '{print $NF}'`

#compute kspace and rpsace data
time $graphene_prog_path/build/data_writer.exe input.ini
#python3 $graphene_prog_path/python/main.py -input input.ini -task prfile

#compute and analyze reciprocal space dynamics
$graphene_prog_path/build/main.exe input.ini
$graphene_prog_path/build/analysis.exe input.ini

#for (( i=1; i<=$Nt; i++ )); do
for i in {$Nt,}; do
    echo "step $i"

    #plot kspace densities
    python3 $graphene_prog_path/python/main.py -input input.ini -tstep $i -task kspace

    #compute and plot rspace densities
    $graphene_prog_path/build/rspace.exe input.ini $i
    python3 $graphene_prog_path/python/main.py -input input.ini -tstep $i -task rspace
done

exit 0
