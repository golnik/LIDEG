import numpy as np
import math

class Graphene:
    def __init__(self,a,s,Z,Nx,Ny):
        self.a = a
        self.s = s
        self.Z = Z

        self.Nx = Nx
        self.Ny = Ny
        
        self.a1 = np.asarray([np.sqrt(3.), 1.])*a/2.
        self.a2 = np.asarray([np.sqrt(3.),-1.])*a/2.
    
        self.A1 = []
        self.A2 = []
        
        zero = np.asarray([0.,0.])
        for ix in range(0,2*Nx+1):
            for iy in range(0,2*Ny+1):
                ixx = ix-Nx
                iyy = iy-Ny

                posA1 = zero + ixx*self.a1 + iyy*self.a2
                posA2 = posA1 - [-self.a/np.sqrt(3.),0.]
                
                self.A1.append(posA1)
                self.A2.append(posA2)
                
    def getA1(self):
        return self.A1
    
    def getA2(self):
        return self.A2
    
    def phi_2pz(self,rx,ry,rz):
        r = np.sqrt(rx**2 + ry**2 + rz**2)
        cosT = rz/r
        return self.Z * r * cosT * np.exp(-0.5*self.Z*r)

    def PhiA(self,rx,ry,rz,kx,ky,A):
        res = 0.
        for posA in A:
            Rx = posA[0]
            Ry = posA[1]
            k_R = kx*Rx + ky*Ry
            res += np.exp(1j*k_R) * self.phi_2pz(rx-Rx,ry-Ry,rz)
        return res
    
    def PhiA1(self,rx,ry,rz,kx,ky):
        return self.PhiA(rx,ry,rz,kx,ky,self.A1)
    
    def PhiA2(self,rx,ry,rz,kx,ky):
        return self.PhiA(rx,ry,rz,kx,ky,self.A2)

    def f(self,kx,ky):
        return np.exp(1j*kx*self.a/np.sqrt(3.)) + 2. * np.exp(-1j*0.5*kx*self.a/np.sqrt(3.)) * np.cos(0.5*ky*self.a)
    
    def phi(self,kx,ky):
        return self.f(kx,ky)/np.abs(self.f(kx,ky))
    
    def psi_p(self,rx,ry,rz,kx,ky):
        return 1./np.sqrt(2.*(1.+self.s*np.abs(self.f(kx,ky)))) \
                *(self.PhiA1(rx,ry,rz,kx,ky) \
                    + np.exp(-1j*self.phi(kx,ky)) * self.PhiA2(rx,ry,rz,kx,ky))