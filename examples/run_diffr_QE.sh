#!/bin/bash

#SBATCH -J graphene
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=90
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=5:00:00
#SBATCH --account=ltamp
#SBATCH --partition=standard

module purge
module add intel/2020.4
module add gnu8/8.3.0
module add openmpi3/3.1.4
module add python/3.9/3.9.10

graphene_prog_path="../"

#input file
input=./inputs/dynamics_input.ini

$graphene_prog_path/build/Extract_QE_WF.exe $input
$graphene_prog_path/build/Transition.exe $input
$graphene_prog_path/build/Diffr_QE.exe $input

python3 $graphene_prog_path/python/comb.py -input $input  #-output plot.pdf


exit 0

