import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

plt.rcParams.update({'figure.figsize': (8,12)})
plt.rcParams.update({'font.size': 16})
plt.rcParams["mathtext.fontset"] = "cm"

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

def plot_tdata(params,tdata,fig_fname):
    fig, ax = plt.subplots(nrows=6, ncols=1, sharex=True)

    tgrid = tdata[:,0]

    ax[0].plot(tgrid,tdata[:,1],label=r"$A_x(t)$")
    ax[0].plot(tgrid,tdata[:,2],label=r"$A_y(t)$")
    ax[0].set_ylabel("Vector\npotential [au]")

    ax[1].plot(tgrid,tdata[:,3],label=r"$E_x(t)$")
    ax[1].plot(tgrid,tdata[:,4],label=r"$E_y(t)$")
    ax[1].set_ylabel("Electric\nfield [au]")

    ax[2].plot(tgrid,tdata[:,5],label=r"$\rho_{vv}(t)$")
    ax[2].plot(tgrid,tdata[:,6],label=r"$\rho_{cc}(t)$")
    ax[2].set_ylabel("Populations")

    ax[3].plot(tgrid,tdata[:,7],label=r"$\mathrm{Re}\{\rho_{cv}(t)\}$")
    ax[3].plot(tgrid,tdata[:,8],label=r"$\mathrm{Im}\{\rho_{cv}(t)\}$")
    ax[3].set_ylabel("Coherences")

    ax[4].plot(tgrid,tdata[:,9] ,label=r"$J_x^{\mathrm{intra}}(t)$")
    ax[4].plot(tgrid,tdata[:,10],label=r"$J_y^{\mathrm{intra}}(t)$")
    ax[4].set_ylabel("Intraband\ncurrent [au]")

    ax[5].plot(tgrid,tdata[:,11],label=r"$J_x^{\mathrm{inter}}(t)$")
    ax[5].plot(tgrid,tdata[:,12],label=r"$J_y^{\mathrm{inter}}(t)$")
    ax[5].set_ylabel("Interband\ncurrent [au]")

    ax[5].set_xlim([params.tmin*au2fs,params.tmax*au2fs])
    ax[5].set_xlabel("Time [fs]")

    for it in range(len(ax)):
        ax[it].legend(loc="upper right")

    fig.align_labels()

    plt.subplots_adjust(left=0.2, bottom=0.1, right=0.95, top=0.97, wspace=0.2, hspace=0.05)

    plt.savefig(fig_fname,dpi=150)
    plt.close()