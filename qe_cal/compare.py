import numpy as np

# Load the two data files
data1 = np.loadtxt('./wf/graphene_wavefunction_1.dat')
data2 = np.loadtxt('./wf/graphene_wavefunction.dat')

# Check if the shapes are the same
if data1.shape != data2.shape:
    print("The files have different shapes.")
else:
    # Compare the content
    if np.array_equal(data1, data2):
        print("The two files are identical.")
    else:
        print("The files are different.")
