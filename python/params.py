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
        self.layers = config['system']['layers']
        self.nlayers = sum(self.layers.count(st) for st in ['A','B','C'])

        self.a = float(config['system']['a'])/au2A
        #self.s = float(config['system']['s'])
        #self.Z = float(config['system']['Z'])

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
        self.dkx = float(config['kgrid']['dkx'])*au2nm
        self.Nkx = int(config['kgrid']['Nkx'])

        self.dky = float(config['kgrid']['dky'])*au2nm
        self.Nky = int(config['kgrid']['Nky'])

        self.Nt = int(config['propagator']['Nt'])

        self.outdir         = config['output']['outdir']
        self.kgfile_fname   = os.path.join(self.outdir,config['output']['kgfile'])
        self.rgfile_fname   = os.path.join(self.outdir,config['output']['rgfile'])
        self.tfile_fname    = os.path.join(self.outdir,config['output']['tfile'])
        self.rhofile_fname  = os.path.join(self.outdir,config['output']['rhofile'])
        self.densfile_fname = os.path.join(self.outdir,config['output']['densfile'])
        self.afile_fname    = os.path.join(self.outdir,config['output']['afile'])
        self.pkfile_fname   = os.path.join(self.outdir,config['output']['pkfile'])
        self.prfile_fname   = os.path.join(self.outdir,config['output']['prfile'])

        #create arrays of atom coordinates
        if os.path.isfile(self.afile_fname):
            with open(self.afile_fname,'r') as file:
                self.atom_coords = []

                Natoms = int(file.readline())
                for ia in range(Natoms):
                    data = file.readline().split()

                    x = float(data[0])
                    y = float(data[1])
                    z = float(data[2])

                    self.atom_coords.append([x,y,z])

        #create grids
        self.tgrid   = np.linspace(self.tmin,self.tmax,self.Nt)

        #we read rgrid from file
        if os.path.isfile(self.rgfile_fname):
            with open(self.rgfile_fname,'r') as file:
                rgrid_data = file.readlines()
            Nxyz = int(rgrid_data[0])

            self.xyzgrid = []

            for iline in range(Nxyz):
                data = rgrid_data[iline+1].split()
                x = float(data[0])
                y = float(data[1])
                z = float(data[2])

                self.xyzgrid.append([x,y,z])

            self.xyzgrid = np.asarray(self.xyzgrid)

        #we read kgrid from file
        if os.path.isfile(self.kgfile_fname):
            with open(self.kgfile_fname,'r') as file:
                kgrid_data = file.readlines()
            Nkxky = int(kgrid_data[0])

            self.kxkygrid = []

            for iline in range(Nkxky):
                data = kgrid_data[iline+1].split()
                kx = float(data[0])
                ky = float(data[1])

                self.kxkygrid.append([kx,ky])

            self.kxkygrid = np.asarray(self.kxkygrid)
