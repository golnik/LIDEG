fname = "input.ini"

isovalues = [-0.0001,0.0001]
#isovalues = [0.001]

col = 2

tsteps = [50]

atom_size = 0.6
bond_size = 0.1
bond_length_trsh = 3.5

###########################################
### Normally no changes after this line ###
###########################################
import sys
import os
import importlib

import numpy as np
import math
from skimage import measure

import bpy
from bpy_extras.object_utils import object_data_add
import bmesh

#path to blender script
script_path = os.path.dirname(bpy.context.space_data.text.filepath)
sys.path.append(script_path)

from UpdateScene import update_scene

#path to graphene python
sys.path.append(os.path.join(script_path,"../python/"))

import params
importlib.reload(params)

def transform_coords(svec,a1,a2):
    #compute angle between a1 and a2 vectors
    ab = np.inner(a1,a2)
    na=np.linalg.norm(a1)
    nb=np.linalg.norm(a2)
    theta=np.arccos(ab/(na*nb))

    #skew-transformation matrix
    S = np.array(((1, np.cos(theta)),(0, np.sin(theta))))

    #rotation matrix
    c,s = np.cos(-theta/2.), np.sin(-theta/2.)
    R = np.array(((c, -s), (s, c)))

    #new xy coordinates are skew and rotated
    vec = R.dot(S.dot([svec[0],svec[1]]))

    #coordinate vectors transformation
    ex = np.asarray([float(params.Nx-1),0.])
    ey = np.asarray([0.,float(params.Ny-1)])

    e_p = ex+ey
    e_m = ex-ey

    s_p = R.dot(S.dot(e_p))
    s_m = R.dot(S.dot(e_m))

    norm_x = np.linalg.norm(a1+a2)/np.linalg.norm(s_p)
    norm_y = np.linalg.norm(a1-a2)/np.linalg.norm(s_m)

    #resulting vector is transformed and normalized
    cvec = [norm_x*vec[0],norm_y*vec[1],svec[2]]

    return cvec

def add_mesh(name, verts, faces, col_name, material):
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(mesh.name, mesh)

    col = bpy.data.collections.get(col_name)
    col.objects.link(obj)

    mesh.from_pydata(verts, [], faces)
    mesh.update()

    for f in mesh.polygons:
        f.use_smooth = True

    mesh.materials.append(material)

    return

def create_dens(col_name,params,tstep,isovalues,materials):
    #load data from prfile
    prdata = np.loadtxt(params.prfile_fname)
    prdata_xyz = prdata[:,0].reshape((params.Nx,params.Ny,params.Nz))

    #load data from rhofile
    rho_t_fname = params.rhofile_fname.replace("%it",'{:06}'.format(tstep))
    rho_data = np.loadtxt(rho_t_fname)

    rho_nocoh_xyz = rho_data[:,0].reshape((params.Nx,params.Ny,params.Nz))
    rho_coh_xyz   = rho_data[:,1].reshape((params.Nx,params.Ny,params.Nz))
    rho_tot_xyz   = rho_data[:,2].reshape((params.Nx,params.Ny,params.Nz))
    
    rhos = [rho_nocoh_xyz-prdata_xyz,    #non coherent
            rho_coh_xyz,                 #coherent
            rho_tot_xyz-prdata_xyz]      #total

    #bravais vectors
    a1 = [params.a/2.*np.sqrt(3.), params.a/2.]
    a2 = [params.a/2.*np.sqrt(3.),-params.a/2.]
    a1 = np.asarray(a1)
    a2 = np.asarray(a2)

    dz = (params.zmax - params.zmin) / (params.Nz - 1) 

    #prepare data for analysis
    Q = rhos[col]
    Qmin = Q.min()
    Qmax = Q.max()

    Niso = len(isovalues)    
    for isv in range(Niso):
        iso = isovalues[isv]
        if iso>Qmin and iso<Qmax:
            #compute vertices and faces for given isovalue
            verts, faces, _, _ = measure.marching_cubes(Q,iso,spacing=(1.,1.,dz),allow_degenerate=False)

            #transform vertices from orthogonal to skew coordinate
            verts_t = []
            for v in verts:
                vt = transform_coords(v,a1,a2)
                verts_t.append(vt)
            verts_t = np.asarray(verts_t)

            verts_all = []
            faces_all = []
            #copy unit cell
            for icx in range(2*params.Nclx+1):
                for icy in range(2*params.Ncly+1):
                    ofset = float(params.Nclx-icx) * a1 + float(params.Ncly-icy) * a2

                    x0 = params.xyzgrid[0,0] - ofset[0]
                    y0 = params.xyzgrid[0,1] - ofset[1]
                    z0 = params.xyzgrid[0,2]

                    Nvertices = len(verts_all)

                    verts_all.extend(verts_t + [x0,y0,z0])
                    faces_all.extend(faces + Nvertices)

            add_mesh(("dens_%03d_%02d" % (tstep,isv)),verts_all,faces_all,col_name,materials[isv])

    return

def create_sphere(name,size,location,col_name,material):
    # Create an empty mesh and the object.
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)
    obj.location = location

    # Add the object into the scene.
    col = bpy.data.collections.get(col_name)
    col.objects.link(obj)

    # Construct the bmesh sphere and assign it to the blender mesh.
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(
        bm,
        u_segments=32,
        v_segments=16,
        radius=size
    )
    bm.to_mesh(mesh)
    bm.free()

    for f in mesh.polygons:
        f.use_smooth = True

    mesh.materials.append(material)
    
    return

def create_bond(x1, y1, z1, x2, y2, z2, r,col_name,material):
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1    
    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    phi = math.atan2(dy, dx) 
    theta = math.acos(dz/dist)

    name = "bond"

    # Create an empty mesh and the object.
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)
    obj.location = (dx/2 + x1, dy/2 + y1, dz/2 + z1)
    obj.rotation_euler[1] = theta 
    obj.rotation_euler[2] = phi 

    # Add the object into the scene.
    col = bpy.data.collections.get(col_name)
    col.objects.link(obj)

    # Select the newly created object
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)

    # Construct the bmesh sphere and assign it to the blender mesh.
    bm = bmesh.new()
    bmesh.ops.create_cone(
        bm,
        cap_ends=True,
        cap_tris=False,
        segments=32,
        radius1=r, 
        radius2=r,
        depth=dist
    )
    bm.to_mesh(mesh)
    bm.free()

    for f in mesh.polygons:
        f.use_smooth = True
        
    mesh.materials.append(material)

    return

def create_graphene(isovalues,params):
    graphene_coll = bpy.data.collections.new("graphene")
    bpy.context.scene.collection.children.link(graphene_coll)

    atoms_coll = bpy.data.collections.new("atoms")
    graphene_coll.children.link(atoms_coll)

    bonds_coll = bpy.data.collections.new("bonds")
    graphene_coll.children.link(bonds_coll)

    dens_coll = bpy.data.collections.new("densities")
    graphene_coll.children.link(dens_coll)
    
    acoords = params.atom_coords
    Natoms = len(acoords)
    
    #apply materials
    atoms_mat = None
    bonds_mat = None
    dens_materials = []
    
    for mat in bpy.data.materials:
        atoms_mat_name     = "atoms_mat"
        bonds_mat_name     = "bonds_mat"
        dens_mat_name_tmp  = "dens_mat_%s"
        
        print(mat.name)
        
        if mat.name == atoms_mat_name:
            atoms_mat = mat
    
        if mat.name == bonds_mat_name:
            bonds_mat = mat

        for iso in range(len(isovalues)):
            dens_mat_name = dens_mat_name_tmp % iso
            if mat.name == dens_mat_name:
                dens_materials.append(mat)
            
    if atoms_mat == None:
        atoms_mat = bpy.data.materials.new(atoms_mat_name)
        atoms_mat.diffuse_color = (1,1,1,1)
    
    if bonds_mat == None:
        bonds_mat = bpy.data.materials.new(bonds_mat_name)
        bonds_mat.diffuse_color = (1,1,1,1)           
    
    if len(dens_materials) == 0:
        max_iso = abs(max(isovalues,key=abs))
        for iso in isovalues:
            dens_mat_name = dens_mat_name_tmp % iso
            dens_mat = bpy.data.materials.new(dens_mat_name)

            color = None
            color_intensity = abs(iso)/max_iso
            print(color_intensity)
            if iso>0:
                color = (color_intensity,0,0,0.5)  #red color for positive iso
            else:
                color = (0,0,color_intensity,0.5)  #blue color for negative iso

            dens_mat.diffuse_color = color
            dens_materials.append(dens_mat)
        
    #create atoms
    print("Create atoms")
    for ia in acoords:
        create_sphere("C",atom_size,(ia[0],ia[1],ia[2]),"atoms",atoms_mat)
    
    #create bonds
    print("Create bonds between atoms")
    for ia in range(Natoms):
        for ja in range(ia+1,Natoms):
            x1 = acoords[ia][0]
            y1 = acoords[ia][1]
            z1 = acoords[ia][2]
            x2 = acoords[ja][0]
            y2 = acoords[ja][1]
            z2 = acoords[ja][2]

            r = np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
            
            if(r<bond_length_trsh):
                create_bond(x1,y1,z1,x2,y2,z2,bond_size,"bonds",bonds_mat)

    #create densities
    print("Create electron densities")

    #if 'tsteps' not in globals():
    
    tsteps = list(range(params.Nt))

    for it in tsteps:
        print("Time step %s" % it)
        iframe = it+1
        create_dens("densities",params,iframe,isovalues,dens_materials)    

    return

if __name__=="__main__":
    params = params.InputParams()

    #read parameters from config file
    if os.path.exists(fname):
        params.analyze_input(fname)
    else:
        print("Input file does not exist!")
    
    #create graphene
    if bpy.data.collections.get("graphene") is None:
        create_graphene(isovalues,params)

    #update scene animation boundaries
    scene = bpy.context.scene 
    scene.frame_end = params.Nt
    
    #update scene to show visible density for specific time step
    update_scene(scene)
    bpy.app.handlers.frame_change_pre.append(update_scene)
