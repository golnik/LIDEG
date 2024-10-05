# Quantum espresso

    ./run_basic.sh

to generate wavefunction

    cd python

    python3 get_all.py

is used to get the C_{n}(k,r) of plane wave saved in coeff folder

# When build exe in hpc

    mkdir build

    cd build

    cmake ../cpp -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc

    make -j

# get density

    cd build

    ./density

# plot 

    cd python

    python3 plot_int_z.py