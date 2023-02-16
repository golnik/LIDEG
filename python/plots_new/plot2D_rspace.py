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

def plot2D_rspace(ax,params,rho_data_xy,cmap,levels):
    if params.rgrid_type == "rectan":
        xgrid = np.linspace(params.xmin,params.xmax,params.Nx)
        ygrid = np.linspace(params.ymin,params.ymax,params.Ny)
        ax.contourf(xgrid*au2A,ygrid*au2A,rho_data_xy,levels=levels,cmap=cmap)
    elif params.rgrid_type == "ucell":
        rho_data_xy = rho_data_xy.flatten(order='C')
        mask = np.isfinite(rho_data_xy)

        xgrid3D = params.xyzgrid[:,0].reshape((params.Nx,params.Ny,params.Nz))
        ygrid3D = params.xyzgrid[:,1].reshape((params.Nx,params.Ny,params.Nz))
        zgrid3D = params.xyzgrid[:,2].reshape((params.Nx,params.Ny,params.Nz))

        xgrid_flat = xgrid3D[:,:,0].flatten(order='C')
        ygrid_flat = ygrid3D[:,:,0].flatten(order='C')

        #bravais vectors
        a1 = [params.a/2.*np.sqrt(3.), params.a/2.]
        a2 = [params.a/2.*np.sqrt(3.),-params.a/2.]
        a1 = np.asarray(a1)
        a2 = np.asarray(a2)

        rho_all = []
        xgrid_all = []
        ygrid_all = []

        for icx in range(-params.Nclx,params.Nclx+1):
            for icy in range(-params.Ncly,params.Ncly+1):
                ofset = float(icx) * a1 + float(icy) * a2
                xgrid = (xgrid_flat[mask] + ofset[0])*au2A
                ygrid = (ygrid_flat[mask] + ofset[1])*au2A

                xgrid_all.extend(xgrid)
                ygrid_all.extend(ygrid)
                rho_all.extend(rho_data_xy[mask])

        cf = ax.tricontourf(xgrid_all,ygrid_all,rho_all,levels=levels,cmap=cmap)

    #acoords = params.atom_coords
    #Natoms = len(acoords)
    #for ia in range(Natoms):
    #    xi = acoords[ia][0]*au2A
    #    yi = acoords[ia][1]*au2A
    #    for ja in range(ia+1,Natoms):
    #        xj = acoords[ja][0]*au2A
    #        yj = acoords[ja][1]*au2A
            
    #        r = np.sqrt((xj-xi)**2 + (yj-yi)**2)

    #        if r<1.5:
    #            ax.plot([xi,xj],[yi,yj],color='black')