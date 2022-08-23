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
    rho_data = rho_data[:].reshape((params.Nx,params.Ny,params.Nz))
    rho_data_xy = np.transpose(np.trapz(rho_data,x=params.zgrid*au2A,axis=2))

    fig, ax = plt.subplots(figsize=(8,8))

    #for ia in params.atom_coords:
    #    x = ia[0]
    #    y = ia[1]
    #    ax.plot(x*au2A,y*au2A,'ro')

    zmin = rho_data_xy.min()
    zmax = rho_data_xy.max()

    #ZZ = max(abs(zmin),abs(zmax))
    ZZ = 6.e-3*0.21
    levels = np.linspace(-ZZ,ZZ,151)

    cf = ax.contourf(params.xgrid*au2A,params.ygrid*au2A,rho_data_xy,levels=levels,cmap=cm.seismic)
    ax.set_box_aspect(1)
    ax.set_xlabel(r"$x [\AA]$",labelpad=10)
    ax.set_ylabel(r"$y [\AA]$",labelpad=5)

    ax.set_xlim([params.xmin*au2A,params.xmax*au2A])
    ax.set_ylim([params.ymin*au2A,params.ymax*au2A])

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(cf, cax=cax, orientation='vertical')

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