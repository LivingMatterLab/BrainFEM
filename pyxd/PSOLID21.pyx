# -*- coding: utf-8 -*-\
# Anisotropic stretch driven growth: Fg = I + (theta-1)*e1*e1, 
# theta_dot = Gs*(stretch-stretchCrit)

from Property cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PSOLID21(Property):
    def __init__(self,material,propProp):
        super(PSOLID21,self).__init__('PSOLID21');
        
        self.Gs           = propProp['Gs'];
        self.stretchCrit  = propProp['stretchCrit'];
        self.material     = material;

    def __str__(self):
        return super(PSOLID21,self).__str__() + "\t Gs = " + str(self.Gs) + "\t stretchCrit = " + str(self.stretchCrit) + "\t material = " + str(self.material.localID)

    cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim):
        cdef int i,j,k,r,l
        cdef double theta, theta0, thRat, stretch, K, dTHdF_frac
        cdef np.ndarray[double, ndim=1] e1
        cdef np.ndarray[double, ndim=2] Fe, Fg, P, dPdTH, dTHdF
        cdef np.ndarray[double, ndim=4] A, Ae
        
        theta0 = Fg0[0,0];
        stretch = np.linalg.norm(F[:,0]);

        if(dt>0):
            theta, K = self.UpdateFg(F,theta0,dt)
        else:
            theta = theta0

        Fg = np.eye(ndim); Fg[0,0] = theta;
        Fginv = np.eye(ndim); Fginv[0,0] = 1./theta;
        Fe = F.dot(Fginv);
        
        Me = self.material.Piola1Stiffness(Fe,ndim);

        e1 = np.zeros(ndim); e1[0] = 1.;
        thRat = (1.-theta)/theta;
        
        P = Me['Pe']
        for i in range(ndim):
            P[i,0]+= thRat*P[i,0];

        if(dt>0):
            Ae = Me['Ae']
            A = np.zeros((ndim,ndim,ndim,ndim))
            for i in range(ndim):
                for j in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim):
                            A[i,j,k,l] = Ae[i,j,k,l] + thRat* (Ae[i,j,k,0]*e1[l] + Ae[i,0,k,l]*e1[j]) \
                                                     + thRat**2 * Ae[i,0,k,0]*e1[j]*e1[l]
            

            # Compute dP_ij/dth = -P_ij/th - Ae_ijrl*Fe_rl / th^2
            dPdTH = np.zeros((ndim,ndim))
            for i in range(ndim):
                for j in range(ndim):
                    dPdTH[i,j] = -P[i,0]*e1[j]/theta
                    for r in range(ndim):
                        dPdTH[i,j]-= Fe[r,0]* (Ae[i,j,r,0]/theta \
                                       + Ae[i,0,r,0]*e1[j]+(1.-theta)/theta**2)

            # Compute dTHdF
            dTHdF_frac = self.Gs*dt/K/theta/stretch
            dTHdF = np.zeros((ndim,ndim))
            for i in range(ndim):
                dTHdF[i,0] = dTHdF_frac*F[i,0]
            
            # Compute A = A + dPdTH_ij*dTHdF_kl
            for i in range(ndim):
                for j in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim):
                            A[i,j,k,l] += dPdTH[i,j]*dTHdF[k,l]
            
            return {'P':P,'A':A,'Fg':Fg}
        else:
            return{'P':P}

    cpdef UpdateFg(self, np.ndarray[double,ndim=2] F,double th0,double dt):
        cdef double theta, stretch, R, K, epsilon
        cdef int it

        theta = th0
        stretch = np.linalg.norm(F[:,0]);

        it = 0;
        epsilon = 1;
        while(epsilon>1e-12 and it<50):
            it+=1;

            R = theta-th0-self.Gs*np.max([stretch/theta-self.stretchCrit,0])*dt;
            K = 1+self.Gs*stretch/theta**2*dt;
            epsilon = np.abs(R/K)

            theta = theta-R/K

        return (theta, K)


