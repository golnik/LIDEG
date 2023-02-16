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

def plot_orb2D(params,it,rho_data,fig_fname):
    rho_data = rho_data.reshape((params.Nx,params.Ny,params.Nz))
    rho_data_xz = np.transpose(np.trapz(rho_data,axis=0))

    fig, ax = plt.subplots(figsize=(8,8))

    xgrid = np.linspace(params.xmin,params.xmax,params.Nx)
    zgrid = np.linspace(params.zmin,params.zmax,params.Nz)
    
    ax.contourf(xgrid*au2A,zgrid*au2A,rho_data_xz,cmap=cm.seismic)

    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.8, top=0.97, wspace=0.2, hspace=0.25)

    plt.savefig(fig_fname,dpi=150)
    plt.close()    

    #print(np.shape(rho_data))

    #rho_data_xz = np.transpose(np.trapz(rho_data,axis=1))