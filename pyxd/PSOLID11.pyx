# -*- coding: utf-8 -*-\
# Anisotropic prescribed growth: Fg = I + (theta-1)*e1*e1
from Property cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PSOLID11(Property):
    def __init__(self,material,propProp):
        super(PSOLID11,self).__init__('PSOLID11');
        
        self.th_rate  = propProp['th_rate'];
        self.material = material;

    def __str__(self):
        return super(PSOLID11,self).__str__() + "\t th_rate = " + str(self.th_rate) + "\t material = " + str(self.material.localID)

    cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim):
        cdef int i,j,k,l
        cdef double theta, thRat
        cdef np.ndarray[double, ndim=1] e1
        cdef np.ndarray[double, ndim=2] Fe, Fg, Fginv, P
        cdef np.ndarray[double, ndim=4] A, Ae
        
        theta = Fg0[0,0]+dt*self.th_rate;
        Fg = np.eye(ndim); Fg[0,0] = theta;
        Fginv = np.eye(ndim); Fginv[0,0] = 1./theta;
        Fe = F.dot(Fginv);
        
        Me = self.material.Piola1Stiffness(Fe,ndim);
        
        e1 = np.zeros(ndim); e1[0] = 1.;
        thRat = (1.-theta)/theta;

        P = Me['Pe']
        for i in range(ndim):
            P[i,0]+= thRat*P[i,0];


        Ae = Me['Ae']
        A = np.zeros((ndim,ndim,ndim,ndim))
        for i in range(ndim):
            for j in range(ndim):
                for k in range(ndim):
                    for l in range(ndim):
                        A[i,j,k,l] = Ae[i,j,k,l] + thRat* (Ae[i,j,k,0]*e1[l] + Ae[i,0,k,l]*e1[j]) \
                                                 + thRat**2 * Ae[i,0,k,0]*e1[j]*e1[l]

        return {'P':P,'A':A,'Fg':Fg}


