import numpy as np
import sys

import plotly.graph_objects as go

def plot_plotly(params,fig_fname):
    pr_data = np.loadtxt(params.prfile_fname)

    dens    = pr_data[:,1]

    X, Y, Z = np.meshgrid(
        np.linspace(params.xmin,params.xmax,params.Nx),
        np.linspace(params.ymin,params.ymax,params.Ny),
        np.linspace(params.zmin,params.zmax,params.Nz)
    )

    data = np.transpose(dens.reshape((params.Nx,params.Ny,params.Nz)),[1,0,2])

    mindata = np.min(data)
    maxdata = np.max(data)

    fig = go.Figure()

    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=data.flatten(),
        isomin=0.01,
        isomax=maxdata,
        opacity=0.1, # max opacity
        opacityscale="uniform",
        surface_count=6,
        colorscale='thermal'
        ))

    xatoms = []
    yatoms = []
    zatoms = []

    for ia in params.atom_coords:
        x = ia[0]
        y = ia[1]
        z = ia[2]

        xatoms.append(x)
        yatoms.append(y)
        zatoms.append(z)

    fig.add_trace(go.Scatter3d(
        mode='markers',
        x=xatoms,
        y=yatoms,
        z=zatoms,
        marker=dict(
            size=12,
            color='black',
            opacity=0.4
        )
        )
    )

    fig.update_layout(
        scene = dict(
            xaxis = dict(nticks=4, range=[params.xmin,params.xmax],),
            yaxis = dict(nticks=4, range=[params.ymin,params.ymax],),
            zaxis = dict(nticks=4, range=[params.zmin,params.zmax],),
        ),
        scene_aspectmode='manual',
        scene_aspectratio=dict(x=1, y=1, z=0.5)
    )

    fig.show()