import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors

plt.rcParams.update({'font.size': 16})
plt.rcParams["mathtext.fontset"] = "cm"

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

def plot_all(params,it,dens_data,rho_data,fig_fname):
    tdata = np.loadtxt(params.tfile_fname)
    tgrid  = tdata[:,0]
    Afield = tdata[:,1]
    Efield = tdata[:,3]
    coh_re = tdata[:,7]

    rho_data = rho_data[:].reshape((params.Nx,params.Ny,params.Nz))
    rho_data_xy = np.transpose(np.trapz(rho_data,x=params.zgrid,axis=2))

    fig = plt.figure(figsize=(9,8))

    gs = fig.add_gridspec(3,2,width_ratios=[1,1], height_ratios=[1,1,3])
    gs.update(wspace=0.5, hspace=0.1)
    axt1 = fig.add_subplot(gs[0,:])
    plt.setp(axt1.get_xticklabels(), visible=False)
    axt2 = fig.add_subplot(gs[1,:],sharex=axt1)
    axl = fig.add_subplot(gs[2,0])
    axr = fig.add_subplot(gs[2,1])

    axt1.plot(tgrid,Afield)
    axt1.set_xlim([params.tmin*au2fs,params.tmax*au2fs])
    axt2.plot(tgrid,coh_re)

    time = params.tgrid[it]*au2fs
    time_str = "Time: %s fs" % ('{:6.2f}'.format(time))

    plt.figtext(0.7,0.93,time_str)

    axt1.axvline(params.tgrid[it]*au2fs,c='r',lw=3)
    axt2.axvline(params.tgrid[it]*au2fs,c='r',lw=3)

    levels = np.linspace(0.,1.,50)
    axl.contourf(params.kx_grid/au2nm,params.ky_grid/au2nm,dens_data,levels=levels,cmap=cm.jet)
    axl.set_box_aspect(1)
    axl.set_xlabel(r"$k_x [\mathrm{nm}^{-1}]$",labelpad=10)
    axl.set_ylabel(r"$k_y [\mathrm{nm}^{-1}]$",labelpad=5)

    rho_m = np.zeros((params.Nx,params.Ny))
    rho_p = np.zeros((params.Nx,params.Ny))
    for ix in range(params.Nx):
        for iy in range(params.Ny):
            rho = rho_data_xy[ix,iy]
            if rho>0:
                rho_p[ix,iy] = rho
            else:
                rho_m[ix,iy] = rho

    #ZZ = max(abs(rho_data_xy.min()),abs(rho_data_xy.max()))
    ZZ = 90.0
    levels_p = np.linspace(0,ZZ,151)
    levels_n = np.linspace(-ZZ,0,151)
    axr.contourf(params.xgrid*au2A,params.ygrid*au2A,rho_p,levels=levels_p,cmap=cm.Reds)
    axr.contourf(params.xgrid*au2A,params.ygrid*au2A,rho_m,levels=levels_n,cmap=cm.Blues)
    axr.set_box_aspect(1)
    axr.set_xlabel(r"$x [\AA]$",labelpad=10)
    axr.set_ylabel(r"$y [\AA]$",labelpad=5)

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
                axr.plot([xi,xj],[yi,yj],color='black')

    #for ia in params.atom_coords:
    #    axr.plot(ia[0]*au2A,ia[1]*au2A,'yo')

    axr.set_xlim([params.xmin*au2A,params.xmax*au2A])
    axr.set_ylim([params.ymin*au2A,params.ymax*au2A])   

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.97, wspace=0.1, hspace=0.1)

    plt.savefig(fig_fname,dpi=150)
    plt.close()