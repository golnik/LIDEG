import argparse
import configparser
import os
import sys
import numpy as np

import params

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

def integrate_multilayer(rho_xyz,params):
    #integrated densities will be stored in this list
    rhos_xy = []

    zgrid_full = np.linspace(params.zmin,params.zmax,params.Nz)

    #case of PH 1layer model (combine 1 layer data to form artificial multilayer)
    if params.model == "hommelhoff":
        #loop over layers
        for il in range(params.nlayers):
            #layer type
            layer_str = params.layers[il]

            sh = 0
            if layer_str == 'A':
                sh = 0
            elif layer_str == 'B':
                sh = 12 #params.a/np.sqrt(3.)
            elif layer_str == 'C':
                sh = 0 #2.*params.a/np.sqrt(3.)
            else:
                raise Exception("Unrecognized layer stacking!")

            rho_xy = np.transpose(np.trapz(rho_xyz,x=zgrid_full,axis=2))

            #shift the 2D array in x and y directions
            rho_xy_sh_x  = np.zeros((params.Nx,params.Ny),dtype=float)
            rho_xy_sh_xy = np.zeros((params.Nx,params.Ny),dtype=float)

            if sh>0:
                #the initial array has a form of ABCDEA where the values at the boundaries are identical
                #we, therefore, create an array  CDEABC according to the specified shift
                rho_xy_sh_x[:sh,:] = rho_xy[params.Nx-sh:,:]
                rho_xy_sh_x[sh:,:] = rho_xy[1:params.Nx-sh+1,:]

                rho_xy_sh_xy[:,:sh] = rho_xy_sh_x[:,params.Ny-sh:]
                rho_xy_sh_xy[:,sh:] = rho_xy_sh_x[:,1:params.Ny-sh+1]
            else:
                rho_xy_sh_xy = rho_xy

            rhos_xy.append(rho_xy_sh_xy)

        #all-layers density contributions are sum over all layers
        rho_xy_all = np.zeros((params.Nx,params.Ny),dtype=float)
        for il in range(params.nlayers):
            rho_xy_all += rhos_xy[il]
            rhos_xy.append(rho_xy_all)

    else: #real multilayer case
        zz_indx = []

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

            zz_indx.append([zindx_min,zindx_max+1])

        zz_indx.append([0,params.Nz])
        nrows = len(zz_indx)

        for iz in range(len(zz_indx)):
            izmin = zz_indx[iz][0]
            izmax = zz_indx[iz][1]

            zgrid = zgrid_full[izmin:izmax]

            rho_xyz_iz = rho_xyz[:,:,izmin:izmax]
            rho_xy = np.transpose(np.trapz(rho_xyz_iz,x=zgrid,axis=2))
            rhos_xy.append(rho_xy)

    return rhos_xy

def plot_tfile(fig,params,tstep,col):
    axs = fig.add_subplot(111)

    tfile_data = np.loadtxt(params.tfile_fname)

    tgrid = tfile_data[:,0]
    data  = tfile_data[:,col]

    axs.plot(tgrid,data)
    axs.set_xlim([params.tmin*au2fs,params.tmax*au2fs])

    time = params.tgrid[tstep-1]*au2fs
    time_str = "Time: %s fs" % ('{:5.2f}'.format(time))

    axs.text(0.75,0.85,time_str,family='monospace',transform=axs.transAxes)

    axs.axvline(time,c='r',lw=3)
    axs.set_xlabel("Time [fs]",labelpad=10)

    plt.subplots_adjust(left=0.12, bottom=0.2, right=0.93, top=0.9)

    return

def plot_rspace(fig,params,rhos,col_names,zmode='auto',layers=True):
    #define number of columns and rows to plot
    ncols = len(col_names)

    nrows = 1
    row_names = ["T"]

    #integrated densities will be stored in this list
    rhos_xy = []

    #integrated density over full z grid
    zgrid_full = np.linspace(params.zmin,params.zmax,params.Nz)
    for icol in range(ncols):
        rho_xyz = rhos[icol]
        rho_xy = np.transpose(np.trapz(rho_xyz,x=zgrid_full,axis=2))
        rhos_xy.append(rho_xy)

    ### multilayer case ###
    if layers == True and params.nlayers > 1:  
        row_names = ['L{}'.format(row+1) for row in range(params.nlayers)] + row_names
        nrows = len(row_names)

        #we clear the integrated densties because in this case they are constructed in a different way
        rhos_xy.clear()

        rhos_xy_ = []
        for ic in range(ncols):
            rhos_xy_ic = integrate_multilayer(rhos[ic],params)
            rhos_xy_.append(rhos_xy_ic)

        for il in range(params.nlayers+1):
            for ic in range(ncols):
                rhos_xy.append(rhos_xy_[ic][il])

    ### plotting ###
    gs = gridspec.GridSpec(nrows, ncols, hspace=0.1, wspace=0.1, left=0.1, bottom=0.2, right=0.95, top=0.9)

    indx = 1 #plot index
    for irow in range(nrows):
        for icol in range(ncols):
            rho_xy = rhos_xy[indx-1]

            Nz = 151
            zmin = rho_xy.min()
            zmax = rho_xy.max()

            if zmode == 'minimax':
                levels = np.linspace(zmin,zmax,Nz)
            elif zmode == 'zabs':
                ZZ = max(abs(zmin),abs(zmax))
                levels = np.linspace(-ZZ,ZZ,Nz)
            else:
                zmin = zmode[0]
                zmax = zmode[1]
                levels = np.linspace(zmin,zmax,Nz)

            ax = fig.add_subplot(gs[irow,icol])

            if irow == 0:
                ax.set_title(col_names[icol])
            if icol == 0:
                ax.set_ylabel(row_names[irow], rotation=0, size='large',labelpad=20)

            plot2D_rspace(ax,params,rho_xy,cm.seismic,levels)

            ax.set_box_aspect(1)
            #ax.set_xlabel(r"$x [\AA]$",labelpad=10)
            #ax.set_ylabel(r"$y [\AA]$",labelpad=5)

            ax.set_xlim([params.xmin*au2A,params.xmax*au2A])
            ax.set_ylim([params.ymin*au2A,params.ymax*au2A])

            #plot_coords(ax,params)

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

def plot_kspace(fig,params,tstep,ist):
    axs = fig.add_subplot(1,1,1)

    dens_t_fname = params.densfile_fname.replace("%it",'{:06}'.format(tstep))

    dens_data = np.loadtxt(dens_t_fname)
    dens_cc = np.transpose(dens_data[:,ist].reshape((params.Nkx,params.Nky)))

    plot2D_kspace(axs,params,dens_cc)#,energies=[1.65,1.65*2])

    axs.set_box_aspect(1)
    axs.set_xlabel(r"$k_x [\mathrm{nm}^{-1}]$",labelpad=10)
    axs.set_ylabel(r"$k_y [\mathrm{nm}^{-1}]$",labelpad=5)

    plt.subplots_adjust(left=0.12, bottom=0.15, right=0.93, top=0.9)

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

    #calculate figure parameters before plotting
    fig_width  = 10
    fig_height = 0
    fig_hratios = np.ones(ntasks)

    for itask in range(ntasks):
        task = tasks[itask]

        if task=="pulse":
            fig_height += 4
            fig_hratios[itask] = 1
        elif task=="kspace":
            fig_height += 5
            fig_hratios[itask] = 2
        elif task=="rspace" or task=="prfile":
            fig_hratios[itask] = 1
            fig_height += 4.0
            if params.nlayers != 1:
                fig_hratios[itask] += params.nlayers
                fig_height += params.nlayers * 2.0
        else:
            raise Exception("Requested task is not available!")

    #plotting
    fig = plt.figure(figsize=(fig_width,fig_height))
    GridSpec = gridspec.GridSpec(ntasks,1,fig,height_ratios=fig_hratios)

    for itask in range(ntasks):
        task = tasks[itask]

        subfig = fig.add_subfigure(GridSpec[itask],frameon=False)

        if task=="pulse":
            plot_tfile(subfig,params,tstep,1)
        elif task=="kspace":
            plot_kspace(subfig,params,tstep,1)
        elif task=="rspace":
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

            #column names
            col_names = ['intra','inter','total']

            plot_rspace(subfig,params,rhos,col_names,zmode='zabs')
        elif task=="prfile":
            rhos = []
            col_names = []

            #load data from prfile
            prdata = np.loadtxt(params.prfile_fname)

            Nst = 2*params.nlayers
            layers = True

            if params.model == "hommelhoff":
                Nst = 2
                layers = False

            for ist in range(Nst):
                prdata_xyz = prdata[:,ist].reshape((params.Nx,params.Ny,params.Nz))
                rhos.append(prdata_xyz)
                col_names.append('S{}'.format(ist+1))

            plot_rspace(subfig,params,rhos,col_names,'minimax',layers)
        else:
            raise Exception("Requested task is not available!")

    plt.savefig(fig_fname,dpi=150,transparent=True)
    plt.close()
    
