import numpy as np
import sys

import plotly.graph_objects as go

def get_data(params,it):
    #load data from prfile
    prdata = np.loadtxt(params.prfile_fname)
    prdata_xyz = prdata[:,0].reshape((params.Nx,params.Ny,params.Nz))

    #load data from rhofile
    rho_t_fname = params.rhofile_fname.replace("%it",'{:06}'.format(it+1))
    rho_data = np.loadtxt(rho_t_fname)

    rho_nocoh_xyz = rho_data[:,0].reshape((params.Nx,params.Ny,params.Nz))
    rho_coh_xyz   = rho_data[:,1].reshape((params.Nx,params.Ny,params.Nz))
    rho_tot_xyz   = rho_data[:,2].reshape((params.Nx,params.Ny,params.Nz))
    
    rhos = [rho_nocoh_xyz-prdata_xyz,
            rho_coh_xyz,
            rho_tot_xyz-prdata_xyz]

    res = rhos[0]
    res = np.transpose(res,[1,0,2])

    return res

def plot3D_rspace(params):
    fr_duration=10  # customize this frame duration according to your data!!!!!
    Nt = params.Nt

    pr_data = np.loadtxt(params.prfile_fname)

    dens    = pr_data[:,1]

    X, Y, Z = np.meshgrid(
        np.linspace(params.xmin,params.xmax,params.Nx),
        np.linspace(params.ymin,params.ymax,params.Ny),
        np.linspace(params.zmin,params.zmax,params.Nz)
    )

    data = np.transpose(dens.reshape((params.Nx,params.Ny,params.Nz)),[1,0,2])

    mindata = -1.e-7
    maxdata = 1.e-7

    fig = go.Figure()

    #plot before the animation starts
    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=get_data(params,0).flatten(),
        isomin=mindata,
        isomax=maxdata,
        opacity=0.1, # max opacity
        opacityscale="extremes",
        surface_count=7,
        colorscale='RdBu'
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

    #create frames
    frames=[
        go.Frame(
            data=[
                go.Volume(
                    x=X.flatten(),
                    y=Y.flatten(),
                    z=Z.flatten(),
                    value=get_data(params,it).flatten(),
                    isomin=mindata,
                    isomax=maxdata,
                    opacity=0.1, # max opacity
                    opacityscale="extremes",
                    surface_count=7,
                    colorscale='RdBu'
                ),
                go.Scatter3d(
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
            ],
            name=f"fr{it}",
            traces=[0,1])
        for it in range(Nt)]

    fig.update(frames=frames)

    def frame_args(duration):
        return {
                "frame": {"duration": duration},
                "mode": "immediate",
                "fromcurrent": True,
                "transition": {"duration": duration, "easing": "linear"},
            }

    sliders = [
                {
                    "pad": {"b": 10, "t": 60},
                    "len": 0.9,
                    "x": 0.1,
                    "y": 0,
                    "steps": [
                        {
                            "args": [[f.name], frame_args(0)],
                            "label": str(k),
                            "method": "animate",
                        }
                        for k, f in enumerate(fig.frames)
                    ],
                }
            ]

    # Layout
    fig.update_layout(
        scene = dict(
            xaxis = dict(nticks=4, range=[params.xmin,params.xmax],),
            yaxis = dict(nticks=4, range=[params.ymin,params.ymax],),
            zaxis = dict(nticks=4, range=[params.zmin,params.zmax],),
        ),
        scene_aspectmode='manual',
        scene_aspectratio=dict(x=1, y=1, z=0.5),
            updatemenus = [
                {
                    "buttons": [
                        {
                            "args": [None, frame_args(50)],
                            "label": "&#9654;", # play symbol
                            "method": "animate",
                        },
                        {
                            "args": [[None], frame_args(0)],
                            "label": "&#9724;", # pause symbol
                            "method": "animate",
                        },
                    ],
                    "direction": "left",
                    "pad": {"r": 10, "t": 70},
                    "type": "buttons",
                    "x": 0.1,
                    "y": 0,
                }
            ],
            sliders=sliders
    )

    fig.show()