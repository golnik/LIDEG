import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

plt.rcParams.update({'figure.figsize': (10,10)})
plt.rcParams.update({'font.size': 16})
plt.rcParams["mathtext.fontset"] = "cm"

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

def plot_debug(params,data,fig_fname):
    fig, axes = plt.subplots(ncols=2)

    data_1 = np.transpose(data[:,0].reshape((params.Nkx,params.Nky)))

    axes[0].contourf(params.kx_grid/au2nm,params.ky_grid/au2nm,data_1,cmap=cm.jet)
    axes[0].set_box_aspect(1)
    axes[0].set_xlabel(r"$k_x [\mathrm{nm}^{-1}]$",labelpad=10)
    axes[0].set_ylabel(r"$k_y [\mathrm{nm}^{-1}]$",labelpad=5)

    data_2 = np.transpose(data[:,1].reshape((params.Nkx,params.Nky)))

    axes[1].contourf(params.kx_grid/au2nm,params.ky_grid/au2nm,data_2,cmap=cm.jet)
    axes[1].set_box_aspect(1)

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.97, wspace=0.2, hspace=0.25)

    plt.savefig(fig_fname,dpi=150)
    plt.close()