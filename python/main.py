import argparse
import configparser
import os
import sys
import numpy as np

import params

import model

from plots.plot_debug import *
from plots.plot_tdata import *
from plots.plot_reciprocal import *
from plots.plot_rspace import *

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

#create figures directory
fig_dir = "figures"
os.makedirs(fig_dir,exist_ok=True)

def plot_tstep(params,it,fig_fname=None):
    if fig_fname==None:
        fig_fname = "fig_%s.png" % ('{:06}'.format(it+1))
        fig_fname = os.path.join(fig_dir,fig_fname)

    print("%s figure is created" % fig_fname)

    #dens_t_fname = params.densfile_fname.replace("%it",'{:06}'.format(it+1))
    #dens_data = np.loadtxt(dens_t_fname)

    #dens_cc = np.transpose(dens_data[:,1].reshape((params.Nkx,params.Nky)))

    #plot_reciprocal(params,it,dens_cc,fig_fname)

    rho_t_fname = params.rhofile_fname.replace("%it",'{:06}'.format(it+1))
    rho_data = np.loadtxt(rho_t_fname)
    rho_data = np.transpose(rho_data[:].reshape((params.Nx,params.Ny)))

    plot_rspace(params,it,rho_data,fig_fname)

    return

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-fname', required=True)
    parser.add_argument('-tstep', required=False)
    parser.add_argument('-out', required=False)

    args = parser.parse_args()
    fname = args.fname

    params = params.InputParams()

    #read parameters from config file
    if os.path.exists(fname):
        params.analyze_input(fname)
    else:
        print("Input file does not exist!")
        sys.exit(1)

    #debug_data = np.loadtxt("debug/debug.out")
    #plot_debug(params,debug_data,"debug/debug.pdf")

    #sys.exit(1)

    #tdata = np.loadtxt(params.tfile_fname)
    #Afield = tdata[:,1]
    #Efield = tdata[:,3]

    #plot time-dependent k-integrated data
    #tfile_fig_name = os.path.join(fig_dir,"tplot.pdf")
    #plot_tdata(params,tdata,tfile_fig_name)

    #plot time-step data
    
    if args.tstep != None:
        it = int(args.tstep)-1
        plot_tstep(params,it,args.out)
    else:
        for it in range(params.Nt):
            plot_tstep(params,it)