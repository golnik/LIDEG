### What should be changed
    In inputs file remember to change,

    K_POINTS {crystal_b}
    (#total k number)

    and remember the time steps!!!!!!!!!!!!!!!!!

### After running the run_kspace.sh

    In the run_para.sh file, please change the t loop from 0 to the number of time steps

    for t in {0..(# of time steps)}; do

    In kpoints.py,

    num_points should be equal to Nkx

### Then run run_para.sh




