# -*- coding: utf-8 -*-\
# Isotropic prescribed growth
from Property cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PSOLID10(Property):
    def __init__(self,material,propProp):
        super(PSOLID10,self).__init__('PSOLID10');
        
        self.th_rate  = propProp['th_rate'];
        self.material = material;

    def __str__(self):
        return super(PSOLID10,self).__str__() + "\t th_rate = " + str(self.th_rate) + "\t material = " + str(self.material.localID)

    def Piola1StiffnessPublic(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim):
        return self.Piola1Stiffness(F,Fg0,dt,ndim)

    cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim):
        cdef double theta
        cdef np.ndarray[double, ndim=2] Fe, Fg, P
        cdef np.ndarray[double, ndim=4] A
        
        theta = Fg0[0,0]+dt*self.th_rate;
        Fe = F/theta;
        Fg = theta*np.eye(ndim);
        
        Me = self.material.Piola1Stiffness(Fe,ndim);
        
        P = Me['Pe']/theta;
        A = Me['Ae']/theta**2;
        
        return {'P':P,'A':A,'Fg':Fg}
