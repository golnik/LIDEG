import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

plt.rcParams.update({'font.size': 16})
plt.rcParams["mathtext.fontset"] = "cm"

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

def plot_rspace(params,it,rho_data,fig_fname):
    rho_data = np.transpose(rho_data[:].reshape((params.Nx,params.Ny)))

    fig, ax = plt.subplots(figsize=(8,8))

    for ia in params.atom_coords:
        x = ia[0]
        y = ia[1]

        ax.plot(x*au2A,y*au2A,'ro')

    zmin = rho_data.min()
    zmax = rho_data.max()

    #ZZ = max(abs(zmin),abs(zmax))
    ZZ = 60.0
    levels = np.linspace(-ZZ,ZZ,50)

    ax.contourf(params.xgrid*au2A,params.ygrid*au2A,rho_data,levels=levels,cmap=cm.seismic)
    ax.set_box_aspect(1)
    ax.set_xlabel(r"$x [\AA]$",labelpad=10)
    ax.set_ylabel(r"$y [\AA]$",labelpad=5)

    ax.set_xlim([params.xmin*au2A,params.xmax*au2A])
    ax.set_ylim([params.ymin*au2A,params.ymax*au2A])

    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.97, wspace=0.2, hspace=0.25)

    plt.savefig(fig_fname,dpi=150)
    plt.close()