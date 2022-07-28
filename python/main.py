import argparse
import configparser
import os
import sys
import numpy as np

import params

import model

from plot_debug import *
from plot_tdata import *
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

    #create figures directory
    fig_dir = "figures"
    os.makedirs(fig_dir,exist_ok=True)

    #debug_data = np.loadtxt("debug/debug.out")
    #plot_debug(params,debug_data,"debug/debug.pdf")

    #sys.exit(1)

    tdata = np.loadtxt(params.tfile_fname)
    Afield = tdata[:,1]
    Efield = tdata[:,3]

    #plot time-dependent k-integrated data
    tfile_fig_name = os.path.join(fig_dir,"tplot.pdf")
    plot_tdata(params,tdata,tfile_fig_name)

    #sys.exit(1)


    #levels = np.linspace(0.,1.,30)

    #dens_00_fname = "output/dens_xy_000001.dat"
    #dens_00_data = np.loadtxt(dens_00_fname)
    #dens_00_data = np.transpose(dens_00_data[:].reshape((Nx,Ny)))

    #gm = model.Graphene(a,s,Z,Nclx,Ncly)
    #A1 = np.asarray(gm.getA1())*au2A
    #A2 = np.asarray(gm.getA2())*au2A

    i = int(params.Nkx/2)
    j = int(params.Nky/2)
    coh_ij = []

    #Nt = 1
    for it in range(params.Nt):
        #it = 49

        print(it)

        dens_t_fname = params.densfile_fname.replace("%it",'{:06}'.format(it+1))
        dens_data  = np.loadtxt(dens_t_fname)
        dens_data = np.transpose(dens_data[:,1].reshape((params.Nkx,params.Nky)))

        #rho_data_re = np.transpose(rho_data[:,2].reshape((params.Nkx,params.Nky)))
        #rho_data_im = np.transpose(rho_data[:,3].reshape((params.Nkx,params.Nky)))

        #rho_data = rho_data_re**2 + rho_data_im**2

        #dens_t_fname = params.densfile_fname.replace("%it",'{:06}'.format(it+1))

        out_t_fname = "fig_%s.png" % ('{:06}'.format(it+1))
        out_t_fname = os.path.join(fig_dir,out_t_fname)
        #out_t_fname = "fig.pdf"

        #coh_ij.append(rho_data[i,j])        

        plot_reciprocal(params,it,Efield,dens_data,out_t_fname)

        #rho_data  = np.loadtxt(rho_t_fname)
        #dens_data = np.loadtxt(dens_t_fname)

        #rho_data = np.transpose(rho_data[:,1].reshape((Nkx,Nky)))
        #dens_data = np.transpose(dens_data[:].reshape((Nx,Ny)))

    sys.exit(1)

    import matplotlib.pyplot as plt
    from matplotlib import cm

    plt.rcParams.update({'figure.figsize': (8,12)})
    plt.rcParams.update({'font.size': 16})
    plt.rcParams["mathtext.fontset"] = "cm"

    fig, ax = plt.subplots()
    ax.plot(tdata[:,0],coh_ij)

    plt.savefig("test.pdf",dpi=150)
    plt.close()