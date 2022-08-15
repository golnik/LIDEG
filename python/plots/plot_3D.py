from mayavi import mlab
import numpy as np
import sys

def plot_3D(params,it,rho_data,fig_fname):
    #x, y, z = np.ogrid[-10:10:20j, -10:10:20j, -10:10:20j]

    x, y, z = np.meshgrid(params.xgrid,params.ygrid,params.zgrid)

    rho_data = rho_data[:].reshape((params.Nx,params.Ny,params.Nz))

    #s = np.sin(x*y*z)/(x*y*z)
    
    #print(x)

    #sys.exit()

    mlab.figure(1, bgcolor=(1, 1, 1), size=(350, 350))
    mlab.clf()

    for ia in params.atom_coords:
        x = ia[0]
        y = ia[1]

        mlab.points3d(x,y,0,
                  scale_factor=1,
                  resolution=20,
                  color=(1, 0, 0),
                  scale_mode='none')

    #s = mlab.mesh(x, y, z)
    
    #mlab.axes(xlabel='x', ylabel='y', zlabel='z',ranges=(0,10000,0,10000,0,22),nb_labels=10)

    #mlab.contour3d(rho_data,extent=[params.xmin,params.xmax,params.ymin,params.ymax,params.zmin,params.zmax])
    #mlab.pipeline.volume(mlab.pipeline.scalar_field(rho_data), figure=fig)
    #mlab.pipeline.volume(mlab.pipeline.scalar_field(rho_data), vmin=0, vmax=0.8)

    #mlab.savefig(filename=fig_fname)

    source = mlab.pipeline.scalar_field(rho_data)
    min = rho_data.min()
    max = rho_data.max()
    vol = mlab.pipeline.volume(source, vmin=min + 0.1 * (max - min),
                                       vmax=min + 0.2 * (max - min))

    mlab.show()