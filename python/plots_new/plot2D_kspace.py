import numpy as np
import sys

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

def plot2D_kspace(ax,params,data,energies=None):
    dens_data = data.flatten(order='F')

    kxkygrid_x = params.kxkygrid[:,0]
    kxkygrid_y = params.kxkygrid[:,1]

    #get position of the Dirac point
    Kx = 2.*np.pi / (np.sqrt(3.)*params.a)
    Ky = 2.*np.pi / (3.*params.a)

    levels = np.linspace(0.,1,100)
    ax.tricontourf((kxkygrid_x-Kx)/au2nm,(kxkygrid_y-Ky)/au2nm,dens_data,levels=levels,cmap=cm.jet)

    #kxgrid = np.linspace(0,1,params.Nkx)
    #kygrid = np.linspace(0,1,params.Nky)
    #ax.contourf(kxgrid/au2nm,kygrid/au2nm,data,levels=levels,cmap=cm.seismic)

    kxmin = -params.dkx
    kxmax = params.dkx
    kymin = -params.dky
    kymax = params.dky

    ax.set_xlim([kxmin/au2nm,kxmax/au2nm])
    ax.set_ylim([kymin/au2nm,kymax/au2nm])

    if energies is not None:
        #load data from pkfile
        pkdata = np.loadtxt(params.pkfile_fname)

        #pkdata = np.loadtxt('/home/u19/ngolubev/calc/graphene/real_field/2.5Vpnm/1layer_PH/output/pkfile.dat')

        gs = 2
        es = 13

        pkdata_E1 = pkdata[:,gs]

        pkdata_E2 = pkdata[:,es]
        dE = pkdata_E2-pkdata_E1

        CS = ax.tricontour((kxkygrid_x-Kx)/au2nm,(kxkygrid_y-Ky)/au2nm,dE*au2eV,levels=energies,linewidths=1.5,colors='red')

        ax.clabel(CS, inline=True, fontsize=10)
