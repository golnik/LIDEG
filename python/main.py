import argparse
import configparser
import os
import sys
import numpy as np

import params

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})
plt.rcParams["mathtext.fontset"] = "cm"

#from plots.plot_debug import *
#from plots.plot_tdata import *
#from plots.plot_reciprocal import *
#from plots.plot_rspace import *
#from plots.plot_2D import *
#from plots.plot_3D import *
#from plots.plot_all import *
#from plots.plot_pyqtgraph import *
#from plots.plot_plotly import *
#from plots.plot_orb2D import *

from plots_new.plot2D_rspace import *
from plots_new.plot2D_kspace import *
#from plots_new.plot3D_rspace import *

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

def find_nearest_indx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def plot_tstep(params,it,fig_fname=None):
    if fig_fname==None:
        fig_fname = "fig_%s.png" % ('{:06}'.format(it+1))
        fig_fname = os.path.join(fig_dir,fig_fname)

    print("%s figure will be created" % fig_fname)

    dens_t_fname = params.densfile_fname.replace("%it",'{:06}'.format(it+1))
    rho_t_fname = params.rhofile_fname.replace("%it",'{:06}'.format(it+1))

    dens_f_exists = os.path.exists(dens_t_fname)
    rho_f_exists  = os.path.exists(rho_t_fname)

    if dens_f_exists:
        dens_data = np.loadtxt(dens_t_fname)
        dens_cc = np.transpose(dens_data[:,1].reshape((params.Nkx,params.Nky)))
        #plot_reciprocal(params,it,dens_cc,fig_fname)

    if rho_f_exists:
        rho_data = np.loadtxt(rho_t_fname)
        #rho_data = np.transpose(rho_data[:].reshape((params.Nx,params.Ny)))

        #plot_rspace(params,it,rho_data,fig_fname)

        #plot_2D(params,it,rho_data,fig_fname)
        #plot_3D(params,it,rho_data,fig_fname)

    #if dens_f_exists and rho_f_exists:
    #    plot_all(params,it,dens_cc,rho_data,fig_fname)

    #plot_plotly(params,fig_fname)

    return

def plot_prfile(params,figname):
    prdata = np.loadtxt(params.prfile_fname)

    fig, axes = plt.subplots(figsize=(8,8),nrows=2,ncols=1)

    for ist in range(2):
        ax=axes[ist]

        prdata_xyz = prdata[:,ist].reshape((params.Nx,params.Ny,params.Nz))
        
        zgrid = np.linspace(params.zmin,params.zmax,params.Nz)
        prdata_xy = np.transpose(np.trapz(prdata_xyz,x=zgrid,axis=2))

        zmin = 0.0
        zmax = 0.22

        levels = np.linspace(zmin,zmax,151)

        plot2D_rspace(ax,params,prdata_xy,cm.jet,levels)

        ax.set_box_aspect(1)
        ax.set_xlabel(r"$x [\AA]$",labelpad=10)
        ax.set_ylabel(r"$y [\AA]$",labelpad=5)

        ax.set_xlim([params.xmin*au2A,params.xmax*au2A])
        ax.set_ylim([params.ymin*au2A,params.ymax*au2A])

    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.8, top=0.97, wspace=0.2, hspace=0.25)

    plt.savefig(figname,dpi=150)
    plt.close()

def plot_rspace(params,tstep,figname,layers=True):
    #load data from prfile
    prdata = np.loadtxt(params.prfile_fname)
    prdata_xyz = prdata[:,0].reshape((params.Nx,params.Ny,params.Nz))

    #load data from rhofile
    rho_t_fname = params.rhofile_fname.replace("%it",'{:06}'.format(tstep))
    rho_data = np.loadtxt(rho_t_fname)

    rho_nocoh_xyz = rho_data[:,0].reshape((params.Nx,params.Ny,params.Nz))
    rho_coh_xyz   = rho_data[:,1].reshape((params.Nx,params.Ny,params.Nz))
    rho_tot_xyz   = rho_data[:,2].reshape((params.Nx,params.Ny,params.Nz))
    
    rhos = [rho_nocoh_xyz-prdata_xyz,
            rho_coh_xyz,
            rho_tot_xyz-prdata_xyz]
 
    zz_indx = []
    
    if layers == True:
        if params.nlayers > 1:
            zvals = np.zeros(params.nlayers+1,dtype=float)
            zvals[0] = params.zmin
            zvals[-1] = params.zmax

            zval = 0.5*params.d
            for il in range(params.nlayers-1):
                zvals[il+1] = zval
                zval += params.d

            for il in range(params.nlayers):
                zindx_min = find_nearest_indx(zgrid_full,zvals[il])
                zindx_max = find_nearest_indx(zgrid_full,zvals[il+1])

                zz_indx.append([zindx_min,zindx_max])

    zz_indx.append([0,params.Nz])
    nrows = len(zz_indx)

    fig, axes = plt.subplots(figsize=(9,nrows*3),nrows=nrows,ncols=3)

    zgrid_full = np.linspace(params.zmin,params.zmax,params.Nz)

    for iz in range(len(zz_indx)):
        izmin = zz_indx[iz][0]
        izmax = zz_indx[iz][1]

        zgrid = zgrid_full[izmin:izmax]

        for ip in range(3):
            rho_xyz = rhos[ip][:,:,izmin:izmax]
            rho_xy = np.transpose(np.trapz(rho_xyz,x=zgrid,axis=2))

            if nrows == 1:
                ax = axes[ip]
            else:
                ax = axes[iz,ip]
            
            zmin = rho_xy.min()
            zmax = rho_xy.max()
            ZZ = max(abs(zmin),abs(zmax))
            #ZZ = 0.01

            levels = np.linspace(-ZZ,ZZ,151)

            plot2D_rspace(ax,params,rho_xy,cm.seismic,levels)

            ax.set_box_aspect(1)
            #ax.set_xlabel(r"$x [\AA]$",labelpad=10)
            #ax.set_ylabel(r"$y [\AA]$",labelpad=5)

            ax.set_xlim([params.xmin*au2A,params.xmax*au2A])
            ax.set_ylim([params.ymin*au2A,params.ymax*au2A])

    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=1, wspace=0.2, hspace=0.25)

    plt.savefig(figname,dpi=150)
    plt.close()    

def plot_kspace(params,tstep,figname,pulse=False):
    Nst = 2*params.nlayers  #get number of states
    ist = Nst-1             #state to plot

    dens_t_fname = params.densfile_fname.replace("%it",'{:06}'.format(tstep))

    dens_data = np.loadtxt(dens_t_fname)
    dens_cc = np.transpose(dens_data[:,ist].reshape((params.Nkx,params.Nky)))

    tfile_data = np.loadtxt(params.tfile_fname)

    fig = plt.figure(figsize=(9,8))

    gs = fig.add_gridspec(2,1,width_ratios=[1],height_ratios=[0.3,1])
    gs.update(wspace=0.5,hspace=0.3)

    axp = fig.add_subplot(gs[0])
    axp.plot(tfile_data[:,0],tfile_data[:,1])
    axp.set_xlim([params.tmin*au2fs,params.tmax*au2fs])

    time = params.tgrid[tstep-1]*au2fs
    time_str = "Time: %s fs" % ('{:6.2f}'.format(time))

    plt.figtext(0.7,0.93,time_str)

    axp.axvline(time,c='r',lw=3)
    axp.set_xlabel("Time [fs]",labelpad=10)

    ax = fig.add_subplot(gs[-1])
    plot2D_kspace(ax,params,dens_cc)

    ax.set_box_aspect(1)
    ax.set_xlabel(r"$k_x [\mathrm{nm}^{-1}]$",labelpad=10)
    ax.set_ylabel(r"$k_y [\mathrm{nm}^{-1}]$",labelpad=5)

    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.97, wspace=0.2, hspace=0.25)

    plt.savefig(figname,dpi=150)
    plt.close()   

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', required=True)
    parser.add_argument('-task', required=True)
    parser.add_argument('-tstep', required=False)
    parser.add_argument('-output', required=False)

    args = parser.parse_args()
    ifname = args.input

    params = params.InputParams()

    #read parameters from config file
    if os.path.exists(ifname):
        params.analyze_input(ifname)
    else:
        print("Input file does not exist!")
        sys.exit(1)

    #create figures directory
    fig_dir = "figures"
    os.makedirs(fig_dir,exist_ok=True)

    fig_fname = args.output
    tstep = 0

    if args.tstep:
        tstep = int(args.tstep)

    if fig_fname==None:
        fig_fname = "fig_%s_%s.png" % (args.task,'{:06}'.format(tstep))
        fig_fname = os.path.join(fig_dir,fig_fname)

    if args.task=="prfile":
        plot_prfile(params,"test.pdf")
    elif args.task=="prfile3D":
        plot_plotly(params,"test.pdf")
    elif args.task=="kspace":
        plot_kspace(params,tstep,fig_fname)
    elif args.task=="rspace":
        plot_rspace(params,tstep,fig_fname)
    elif args.task=="rspace3D":
        plot3D_rspace(params)
    else:
        raise Exception("Requested task is not available!")
    