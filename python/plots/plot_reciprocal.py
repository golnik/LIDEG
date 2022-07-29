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

def plot_reciprocal(params,it,dens_data,fig_fname):
    fig, ax = plt.subplots(figsize=(8,8))

    levels = np.linspace(0.,1.,100)

    ax.contourf(params.kx_grid/au2nm,params.ky_grid/au2nm,dens_data,levels=levels,cmap=cm.jet)
    ax.set_box_aspect(1)
    ax.set_xlabel(r"$k_x [\mathrm{nm}^{-1}]$",labelpad=10)
    ax.set_ylabel(r"$k_y [\mathrm{nm}^{-1}]$",labelpad=5)

    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.97, wspace=0.2, hspace=0.25)

    plt.savefig(fig_fname,dpi=150)
    plt.close()

def plot_reciprocal_with_pulse(params,it,Efield_data,rho_data,fig_fname):
    fig, axes = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [1,3]})

    axes[0].plot(params.tgrid*au2fs,Efield_data*au2Vnm)
    axes[0].set_xlim([params.tmin*au2fs,params.tmax*au2fs])
    axes[0].set_xlabel("Time [fs]",labelpad=10)
    axes[0].set_ylabel("E [V/nm]")

    time = params.tgrid[it]*au2fs
    time_str = "Time: %s fs" % ('{:6.2f}'.format(time))

    plt.figtext(0.7,0.93,time_str)

    axes[0].axvline(params.tgrid[it]*au2fs,c='r',lw=3)

    levels = np.linspace(0.,1.,100)
    #levels = np.linspace(-0.5,0.5,100)

    axes[1].contourf(params.kx_grid/au2nm,params.ky_grid/au2nm,rho_data,levels=levels,cmap=cm.jet)
    axes[1].set_box_aspect(1)
    axes[1].set_xlabel(r"$k_x [\mathrm{nm}^{-1}]$",labelpad=10)
    axes[1].set_ylabel(r"$k_y [\mathrm{nm}^{-1}]$",labelpad=5)

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.97, wspace=0.2, hspace=0.25)

    plt.savefig(fig_fname,dpi=150)
    plt.close()