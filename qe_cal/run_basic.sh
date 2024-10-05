#!/bin/bash

module add qe/7.0

# module add compiler/2024.0.2
# module add mkl/2024.0
# module add mpi/2021.11
# module add hdf5/intel2024.0.2

# export PATH=$PATH:"/home/mingruiyuan/software/test_qe/q-e-qe-7.2/build/bin"


rm ./coeff ./wf ./figures ./tmp ./outputs -rf

mkdir ./coeff ./wf ./figures ./outputs

mpirun -np 8 pw.x -inp inputs/graphene_scf.in > outputs/graphene_scf.out

#mpirun -np 8 pw.x -in inputs/graphene_nscf.in > outputs/graphene_nscf.out

mpirun -np 8 dos.x -inp inputs/graphene_dos.in > outputs/graphene_dos.out

mpirun -np 8 pw.x -inp inputs/graphene_bands.in > outputs/graphene_bands.out

mpirun -np 8 bands.x -inp inputs/graphene_bands_pp.in > outputs/graphene_bands_pp.out

exit 0
