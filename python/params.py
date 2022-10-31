import configparser
import numpy as np
import sys
import os

au2nm  = 0.052917829614246
au2A   = au2nm*10
au2eV  = 27.211396641308
au2Vnm = 5.14220826*10**2
au2fs  = 0.02418884254

class InputParams:
    def __init__(self):
        self.config = configparser.ConfigParser(interpolation=None)

    def analyze_input(self,fname):
        config = self.config
        config.read(fname)

        #read model parameters
        self.a = float(config['model']['a'])/au2A
        self.s = float(config['model']['s'])
        self.Z = float(config['model']['Z'])

        #read tgrid
        self.tmin = float(config['propagator']['tmin'])/au2fs
        self.tmax = float(config['propagator']['tmax'])/au2fs
        self.Nt = int(config['propagator']['Nt'])

        #read rgrid
        self.rgrid_type = config['rgrid']['type']
        self.Nclx = int(config['rgrid']['Nclx'])
        self.Ncly = int(config['rgrid']['Ncly'])

        self.xmin = float(config['rgrid']['xmin'])/au2A
        self.xmax = float(config['rgrid']['xmax'])/au2A
        self.Nx = int(config['rgrid']['Nx'])

        self.ymin = float(config['rgrid']['ymin'])/au2A
        self.ymax = float(config['rgrid']['ymax'])/au2A
        self.Ny = int(config['rgrid']['Ny'])

        self.zmin = float(config['rgrid']['zmin'])/au2A
        self.zmax = float(config['rgrid']['zmax'])/au2A       
        self.Nz = int(config['rgrid']['Nz'])

        #read kgrid
        #self.dkx = float(config['kgrid']['dkx'])*au2nm
        self.Nkx = int(config['kgrid']['Nkx'])

        #self.dky = float(config['kgrid']['dky'])*au2nm
        self.Nky = int(config['kgrid']['Nky'])

        self.Nt = int(config['propagator']['Nt'])

        self.outdir         = config['output']['outdir']
        self.kgfile_fname   = os.path.join(self.outdir,config['output']['kgfile'])
        self.rgfile_fname   = os.path.join(self.outdir,config['output']['rgfile'])
        self.tfile_fname    = os.path.join(self.outdir,config['output']['tfile'])
        self.rhofile_fname  = os.path.join(self.outdir,config['output']['rhofile'])
        self.densfile_fname = os.path.join(self.outdir,config['output']['densfile'])
        self.afile_fname    = os.path.join(self.outdir,config['output']['afile'])

        #create arrays of atom coordinates
        if os.path.isfile(self.afile_fname):
            with open(self.afile_fname,'r') as file:
                self.atom_coords = []

                Natoms = int(file.readline())
                for ia in range(Natoms):
                    data = file.readline().split()

                    #A1 atom
                    x = float(data[0])
                    y = float(data[1])
                    self.atom_coords.append([x,y])

                    #A2 atom
                    x = float(data[2])
                    y = float(data[3])
                    self.atom_coords.append([x,y])

        #create grids
        self.tgrid   = np.linspace(self.tmin,self.tmax,self.Nt)

        #we read rgrid from file
        if os.path.isfile(self.rgfile_fname):
            with open(self.rgfile_fname,'r') as file:
                rgrid_data = file.readlines()
            Nxy = int(rgrid_data[0])

            self.xgrid = []
            self.ygrid = []

            for iline in range(Nxy):
                data = rgrid_data[iline+1].split()
                x = float(data[0])
                y = float(data[1])

                self.xgrid.append(x)
                self.ygrid.append(y)

            self.xgrid = np.asarray(self.xgrid)
            self.ygrid = np.asarray(self.ygrid)

            #self.xgrid = np.linspace(self.xmin,self.xmax,self.Nx)
            #self.ygrid = np.linspace(self.ymin,self.ymax,self.Ny)
            self.zgrid = np.linspace(self.zmin,self.zmax,self.Nz)

        #we read kgrid from file
        if os.path.isfile(self.kgfile_fname):
            with open(self.kgfile_fname,'r') as file:
                kgrid_data = file.readlines()
            Nkxky = int(kgrid_data[0])

            self.kx_grid = []
            self.ky_grid = []

            for iline in range(Nkxky):
                data = kgrid_data[iline+1].split()
                kx = float(data[0])
                ky = float(data[1])

                self.kx_grid.append(kx)
                self.ky_grid.append(ky)

            self.kx_grid = np.asarray(self.kx_grid)
            self.ky_grid = np.asarray(self.ky_grid)            

            #self.kx_grid = np.float_(kgrid_data[1].split())
            #self.ky_grid = np.float_(kgrid_data[3].split())