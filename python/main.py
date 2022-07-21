import argparse
import configparser
import os
import sys
import numpy as np

import params

import model

from plot_reciprocal import *

#import matplotlib.pyplot as plt
#from matplotlib import cm

#plt.rcParams.update({'figure.figsize': (12,8)})
#plt.rcParams.update({'font.size': 16})
#plt.rcParams["mathtext.fontset"] = "cm"

#plt.rcParams['text.usetex'] = True
#plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-fname')

    args = parser.parse_args()
    fname = args.fname

    params = params.InputParams()

    #read parameters from config file
    if os.path.exists(fname):
        params.analyze_input(fname)
    else:
        print("Input file does not exist!")
        sys.exit(1)

    tdata = np.loadtxt(params.tfile_fname)
    Afield = tdata[:,1]
    Efield = tdata[:,3]

    levels = np.linspace(0.,1.,30)

    #dens_00_fname = "output/dens_xy_000001.dat"
    #dens_00_data = np.loadtxt(dens_00_fname)
    #dens_00_data = np.transpose(dens_00_data[:].reshape((Nx,Ny)))

    #gm = model.Graphene(a,s,Z,Nclx,Ncly)
    #A1 = np.asarray(gm.getA1())*au2A
    #A2 = np.asarray(gm.getA2())*au2A
    
    #Nt = 1
    for it in range(params.Nt):
        #it = 150

        print(it)

        rho_t_fname = params.rhofile_fname.replace("%it",'{:06}'.format(it+1))
        rho_data  = np.loadtxt(rho_t_fname)
        rho_data = np.transpose(rho_data[:,1].reshape((params.Nkx,params.Nky)))

        dens_t_fname = params.densfile_fname.replace("%it",'{:06}'.format(it+1))

        out_t_fname = "output/fig_%s.png" % ('{:06}'.format(it+1))
        #out_t_fname = "fig.pdf"

        plot_reciprocal(params,it,Efield,rho_data,out_t_fname)

        #rho_data  = np.loadtxt(rho_t_fname)
        #dens_data = np.loadtxt(dens_t_fname)

        #rho_data = np.transpose(rho_data[:,1].reshape((Nkx,Nky)))
        #dens_data = np.transpose(dens_data[:].reshape((Nx,Ny)))

