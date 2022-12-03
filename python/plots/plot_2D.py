import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams.update({'font.size': 16})
plt.rcParams["mathtext.fontset"] = "cm"

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

def plot_2D(params,it,rho_data,fig_fname):
    if params.Nz == 0 or params.Nz == 1:
        rho_data_vv    = rho_data[:,0].reshape((params.Nx,params.Ny,1))
        rho_data_cc    = rho_data[:,1].reshape((params.Nx,params.Ny,1))
        rho_data_nocoh = rho_data[:,2].reshape((params.Nx,params.Ny,1))
        rho_data_coh   = rho_data[:,3].reshape((params.Nx,params.Ny,1))
        rho_data_total = rho_data[:,4].reshape((params.Nx,params.Ny,1))

        rho_data_xy = rho_data_total - rho_data_vv
        #rho_data_xy = rho_data_coh
    #else:
    #    rho_data = rho_data[:,col].reshape((params.Nx,params.Ny,params.Nz))
    #    rho_data_xy = np.transpose(np.trapz(rho_data,x=params.zgrid*au2A,axis=2))

    fig, ax = plt.subplots(figsize=(8,8))

    #for ia in params.atom_coords:
    #    x = ia[0]
    #    y = ia[1]
    #    ax.plot(x*au2A,y*au2A,'ro')

    zmin = rho_data_xy.min()
    zmax = rho_data_xy.max()

    ZZ = max(abs(zmin),abs(zmax))
    #ZZ = 0.001
    #print(ZZ)
    levels = np.linspace(-ZZ,ZZ,151)

    if params.rgrid_type == "regular":
        xgrid = np.linspace(params.xmin,params.xmax,params.Nx)
        ygrid = np.linspace(params.ymin,params.ymax,params.Ny)
        ax.contourf(xgrid*au2A,ygrid*au2A,rho_data_xy,levels=levels,cmap=cm.seismic)
    elif params.rgrid_type == "ucell":
        rho_data_xy = rho_data_xy.flatten(order='C')
        mask = np.isfinite(rho_data_xy)

        #bravais vectors
        a1 = [params.a/2.*np.sqrt(3.), params.a/2.]
        a2 = [params.a/2.*np.sqrt(3.),-params.a/2.]
        a1 = np.asarray(a1)
        a2 = np.asarray(a2)

        rho_all = []
        xgrid_all = []
        ygrid_all = []

        #define the required number of unit cells
        Nclx = [2,2]
        Ncly = [2,1]

        for icx in range(-Nclx[0],Nclx[1]+1):
            for icy in range(-Ncly[0],Ncly[1]+1):
                ofset = float(icx) * a1 + float(icy) * a2
                xgrid = (params.xgrid[mask] + ofset[0])*au2A
                ygrid = (params.ygrid[mask] + ofset[1])*au2A

                xgrid_all.extend(xgrid)
                ygrid_all.extend(ygrid)
                rho_all.extend(rho_data_xy[mask])

        cf = ax.tricontourf(xgrid_all,ygrid_all,rho_all,levels=levels,cmap=cm.seismic)

    ax.set_box_aspect(1)
    ax.set_xlabel(r"$x [\AA]$",labelpad=10)
    ax.set_ylabel(r"$y [\AA]$",labelpad=5)

    ax.set_xlim([params.xmin*au2A,params.xmax*au2A])
    ax.set_ylim([params.ymin*au2A,params.ymax*au2A])

    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes('right', size='5%', pad=0.05)
    #fig.colorbar(cf, cax=cax, orientation='vertical')

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

    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.8, top=0.97, wspace=0.2, hspace=0.25)

    plt.savefig(fig_fname,dpi=150)
    plt.close()