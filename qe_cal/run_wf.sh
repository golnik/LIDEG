#!/bin/bash

module add qe/7.0

mkdir -p ./wf

mpirun -np 8 pp.x < inputs/graphene_2D_wf.in > outputs/graphene_2D_wf.out

#mpirun -np 4 pp.x < inputs/pp_wavefunction_2D.in > outputs/pp_wavefunction_2D.out