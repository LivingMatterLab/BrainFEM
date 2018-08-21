# -*- coding: utf-8 -*-\
# Isotropic stretch driven growth
from Property cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PSOLID20(Property):
    def __init__(self,material,propProp):
        super(PSOLID20,self).__init__('PSOLID20');
        
        self.Gs     = propProp['Gs'];
        self.Jcrit  = propProp['Jcrit'];
        self.material = material;

    def __str__(self):
        return super(PSOLID20,self).__str__() + "\t Gs = " + str(self.Gs) + "\t Jcrit = " + str(self.Jcrit) + "\t material = " + str(self.material.localID)

    cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim):
        cdef int i,j,k,r,l
        cdef double theta, theta0, K
        cdef np.ndarray[double, ndim=2] Finv, Fe, Fg, P, dPdTH, dTHdF
        cdef np.ndarray[double, ndim=4] A, Ae
        
        theta0 = Fg0[0,0];

        if(dt>0):
            theta, K = self.UpdateFg(F,theta0,dt,ndim)
        else:
            theta = theta0

        Fe = F/theta;
        Fg = theta*np.eye(ndim);
        
        Me = self.material.Piola1Stiffness(Fe,ndim);
        
        P = Me['Pe']/theta;

        if(dt>0):
            Ae = Me['Ae'];

            # Compute dP_ij/dth = -P_ij/th - Ae_ijrl*Fe_rl / th^2
            dPdTH = -P/theta;
            for i in range(ndim):
                for j in range(ndim):
                    for r in range(ndim):
                        for l in range(ndim):
                            dPdTH[i,j] -= Ae[i,j,r,l]*Fe[r,l]/theta**2


            # Compute dTHdF
            Finv = np.linalg.inv(F)
            dTHdF = self.Gs*dt/K*np.linalg.det(F)*Finv.T/theta**ndim;
            #dTHdF = self.Gs*dt/K*np.linalg.det(F)*Finv.T/2/theta**3;
            
            print dPdTH
            print dTHdF


            # Compute A = dPdF_ijkl + dPdTH_ij*dTHdF_kl
            A = np.zeros((ndim,ndim,ndim,ndim));
            for i in range(ndim):
                for j in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim):
                            A[i,j,k,l] = Ae[i,j,k,l]/theta**2 + dPdTH[i,j]*dTHdF[k,l]
            
            return {'P':P,'A':A,'Fg':Fg}
        else:
            return{'P':P}

    cpdef UpdateFg(self, np.ndarray[double,ndim=2] F,double th0,double dt, int ndim):
        cdef double theta, J, R, K, epsilon
        cdef int it

        theta = th0
        J = np.linalg.det(F)

        it = 0;
        epsilon = 1;
        while(epsilon>1e-12 and it<50):
            it+=1;

            R = theta-th0-self.Gs*np.max([J/theta**ndim-self.Jcrit,0])*dt;
            K = 1+ndim*self.Gs*J/theta**(ndim+1)*dt;
            epsilon = np.abs(R/K)

            theta = theta-R/K

        return (theta, K)


