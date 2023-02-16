import numpy as np

import pyqtgraph.examples

import pyqtgraph as pg
import pyqtgraph.opengl as gl

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

def plot_pyqtgraph(params,fig_fname):
    #pyqtgraph.examples.run()

    app = pg.mkQApp("GLIsosurface Example")
    w = gl.GLViewWidget()
    w.show()
    w.setWindowTitle('pyqtgraph example: GLIsosurface')

    w.setCameraPosition(distance=40)

    g = gl.GLGridItem()
    g.scale(2,2,1)
    w.addItem(g)

    md = gl.MeshData.sphere(rows=10, cols=20, radius=0.5)

    for atom in params.atom_coords:
        m1 = gl.GLMeshItem(
            meshdata=md,
            smooth=True,
            color=(1,0,0,0.5),
            shader="balloon",
            glOptions="additive",
        )
        m1.translate(atom[0],atom[1],atom[2])

        w.addItem(m1)

    pr_data = np.loadtxt(params.prfile_fname)
    
    xcoords = pr_data[:,0]
    ycoords = pr_data[:,1]
    zcoords = pr_data[:,2]
    dens    = pr_data[:,3]

    data = dens.reshape((params.Nx,params.Ny,params.Nz))

    rho_data_xy = np.transpose(np.trapz(data,axis=2))

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
        Nclx = [0,0]
        Ncly = [0,0]

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


    #d2 = np.empty(data.shape + (4,), dtype=np.ubyte)

    #positive = np.log(np.clip(data, 0, data.max())**2)
    #negative = np.log(np.clip(-data, 0, -data.min())**2)

    #d2 = np.empty(data.shape + (4,), dtype=np.ubyte)
    #d2[..., 0] = positive * (255./positive.max())
    #d2[..., 1] = negative * (255./negative.max())
    #d2[..., 2] = d2[...,1]
    #d2[..., 3] = d2[..., 0]*0.3 + d2[..., 1]*0.3
    #d2[..., 3] = (d2[..., 3].astype(float) / 255.) **2 * 255

    #d2[:, 0, 0] = [255,0,0,100]
    #d2[0, :, 0] = [0,255,0,100]
    #d2[0, 0, :] = [0,0,255,100]

    #v = gl.GLVolumeItem(d2)
    #v.translate(-50,-50,-100)
    #w.addItem(v)

    #verts, faces = pg.isosurface(data,0.001)

    #md = gl.MeshData(vertexes=verts, faces=faces)

    #colors = np.ones((md.faceCount(), 4), dtype=float)
    #colors[:,3] = 0.2
    #colors[:,2] = np.linspace(0, 1, colors.shape[0])
    #md.setFaceColors(colors)
    #m1 = gl.GLMeshItem(meshdata=md, smooth=False, shader='balloon')
    #m1.setGLOptions('additive')

    #w.addItem(m1)
    #m1.translate(-25, -25, -20)    
    
    #dens_cc = np.transpose(dens_data[:,3].reshape((params.Nkx,params.Nky)))

    #pg.exec()