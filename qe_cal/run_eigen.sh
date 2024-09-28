#!/bin/bash

module add qe/7.0

# module add compiler/2024.0.2
# module add mkl/2024.0
# module add mpi/2021.11
# module add hdf5/intel2024.0.2

# export PATH=$PATH:"/home/mingruiyuan/software/test_qe/q-e-qe-7.2/build/bin"

mkdir -p ./eigens

mpirun -np 8 projwfc.x < ./inputs/projwfc.in > ./outputs/projwfc.out

