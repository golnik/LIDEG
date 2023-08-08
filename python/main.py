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

def plot_tfile(fig,params,tstep,col):
    axs = fig.add_subplot(111)

    tfile_data = np.loadtxt(params.tfile_fname)

    tgrid = tfile_data[:,0]
    data  = tfile_data[:,col]

    axs.plot(tgrid,data)
    axs.set_xlim([params.tmin*au2fs,params.tmax*au2fs])

    time = params.tgrid[tstep-1]*au2fs
    time_str = "Time: %s fs" % ('{:6.2f}'.format(time))

    axs.text(0.7,0.9,time_str,transform=plt.gcf().transFigure)

    axs.axvline(time,c='r',lw=3)
    axs.set_xlabel("Time [fs]",labelpad=10)

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

def plot_rspace(fig,params,tstep,layers=False):
    #load data from prfile
    prdata = np.loadtxt(params.prfile_fname)
    prdata_xyz = prdata[:,0].reshape((params.Nx,params.Ny,params.Nz))

    #load data from rhofile
    rho_t_fname = params.rhofile_fname.replace("%it",'{:06}'.format(tstep))
    rho_data = np.loadtxt(rho_t_fname)

    #reshape densities to 3D arrays
    rho_nocoh_xyz = rho_data[:,0].reshape((params.Nx,params.Ny,params.Nz))
    rho_coh_xyz   = rho_data[:,1].reshape((params.Nx,params.Ny,params.Nz))
    rho_tot_xyz   = rho_data[:,2].reshape((params.Nx,params.Ny,params.Nz))
    
    #densities for plotting
    rhos = [rho_nocoh_xyz-prdata_xyz,
            rho_coh_xyz,
            rho_tot_xyz-prdata_xyz]
 
    #define number of columns and rows to plot
    ncols = 3
    col_names = ['intra','inter','total']

    nrows = 1
    row_names = []
    if layers == True:
        row_names = ['L{}'.format(row+1) for row in range(params.nlayers)]
    row_names.append("T")

    #create required z grids for multilayer plotting
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

    #overall grid
    zz_indx.append([0,params.Nz])

    zgrid_full = np.linspace(params.zmin,params.zmax,params.Nz)

    indx = 1 #plot index
    for iz in range(len(zz_indx)):
        izmin = zz_indx[iz][0]
        izmax = zz_indx[iz][1]

        zgrid = zgrid_full[izmin:izmax]

        #plot nocoh, coh, and total densities
        for ip in range(ncols):
            rho_xyz = rhos[ip][:,:,izmin:izmax]
            rho_xy = np.transpose(np.trapz(rho_xyz,x=zgrid,axis=2))

            zmin = rho_xy.min()
            zmax = rho_xy.max()
            ZZ = max(abs(zmin),abs(zmax))
            #ZZ = 0.01

            levels = np.linspace(-ZZ,ZZ,151)

            ax = fig.add_subplot(nrows,ncols,indx)

            if iz == 0:
                ax.set_title(col_names[ip])
            if ip == 0:
                ax.set_ylabel(row_names[iz], rotation=0, size='large',labelpad=20)

            plot2D_rspace(ax,params,rho_xy,cm.seismic,levels)

            ax.set_box_aspect(1)
            #ax.set_xlabel(r"$x [\AA]$",labelpad=10)
            #ax.set_ylabel(r"$y [\AA]$",labelpad=5)

            ax.set_xlim([params.xmin*au2A,params.xmax*au2A])
            ax.set_ylim([params.ymin*au2A,params.ymax*au2A])

            plot_coords(ax,params)

            indx += 1

    return

def plot_coords(ax,params):
    acoords = params.atom_coords
    Natoms = len(acoords)
    for ia in range(Natoms):
        xi = acoords[ia][0]*au2A
        yi = acoords[ia][1]*au2A
        for ja in range(ia+1,Natoms):
            xj = acoords[ja][0]*au2A
            yj = acoords[ja][1]*au2A
            
            r = np.sqrt((xj-xi)**2 + (yj-yi)**2)

            if r<1.5:
                ax.plot([xi,xj],[yi,yj],color='black')
    return

def plot_kspace(fig,params,tstep):
    axs = fig.add_subplot(1,1,1)

    Nst = 2*params.nlayers  #get number of states
    ist = Nst-1             #state to plot

    dens_t_fname = params.densfile_fname.replace("%it",'{:06}'.format(tstep))

    dens_data = np.loadtxt(dens_t_fname)
    dens_cc = np.transpose(dens_data[:,ist].reshape((params.Nkx,params.Nky)))

    plot2D_kspace(axs,params,dens_cc)#,energies=[1.65,1.65*2])

    axs.set_box_aspect(1)
    axs.set_xlabel(r"$k_x [\mathrm{nm}^{-1}]$",labelpad=10)
    axs.set_ylabel(r"$k_y [\mathrm{nm}^{-1}]$",labelpad=5)

    return

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-input', required=True)
    parser.add_argument('-tasks', required=True)
    parser.add_argument('-tstep', required=True)
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

    tasks_str = args.tasks
    tasks = args.tasks.split(",")
    ntasks = len(tasks)

    if ntasks == 0:
        print("Tasks are not specified!")
        sys.exit(1)

    #create figures directory
    fig_dir = "figures"
    os.makedirs(fig_dir,exist_ok=True)

    fig_fname = args.output
    tstep = int(args.tstep)

    if fig_fname==None:
        fig_fname = "fig_%s_%s.png" % (tasks_str,'{:06}'.format(tstep))
        fig_fname = os.path.join(fig_dir,fig_fname)

    #create figure and required axes
    width  = 10
    height = 4.5 * ntasks

    fig = plt.figure(figsize=(width,height))
    subfigs = fig.subfigures(ntasks, 1)

    for itask in range(ntasks):
        task = tasks[itask]

        if ntasks == 1:
            subfig = subfigs
        else:
            subfig = subfigs[itask]

        if task=="pulse":
            plot_tfile(subfig,params,tstep,1)
        elif task=="kspace":
            plot_kspace(subfig,params,tstep)
        elif task=="rspace":
            plot_rspace(subfig,params,tstep)
        else:
            raise Exception("Requested task is not available!")

    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.97, wspace=0.2, hspace=0.5)

    plt.savefig(fig_fname,dpi=150)
    plt.close()
    
