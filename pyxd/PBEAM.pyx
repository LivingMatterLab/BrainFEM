# -*- coding: utf-8 -*-\
# Bar element, assumed to be incompressible
from Property cimport *
from DataContainer cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PBEAM(Property):
    def __init__(self,material,propProp):
        super(PBEAM,self).__init__('PBEAM');
        
        self.area  = propProp['area'];
        self.I11   = propProp['I11'];
        self.I22   = propProp['I22'];
        self.I12   = propProp['I12'];
        self.material = material;

    def __str__(self):
        return super(PBEAM,self).__str__() + "\t Area = " + str(self.area) + "\t [I11,I22,I12] = [" + str(self.I11) +","+ str(self.I22) +","+ str(self.I12) +"]\t material = " + str(self.material.localID)

    cpdef Piola1Stiffness(self,np.ndarray[double,ndim=1] F, np.ndarray[double,ndim=1] curvature,double dt,int ndim):
        cdef np.ndarray[double,ndim=1] P,M
        cdef np.ndarray[double,ndim=2] dPdF, dMdK

        # Force and moment vector
        P = np.array([self.material.E *self.area*(F[0]),\
                      self.material.mu*self.area*F[1],\
                      self.material.mu*self.area*F[2]])
        
        M = np.array([self.material.mu*curvature[0]*(self.I11+self.I22),\
                      self.material.E *curvature[1]*self.I11,\
                      self.material.E *curvature[2]*self.I22])

        # Material tangents
        dPdF = np.diag(np.array([self.material.E *self.area, \
                                 self.material.mu*self.area, \
                                 self.material.mu*self.area]))
        
        dMdK = np.diag(np.array([self.material.mu*(self.I11+self.I22),\
                                 self.material.E *self.I11,\
                                 self.material.E *self.I22]))

        return {'P':P,'M':M,'dPdF':dPdF,'dMdK':dMdK}

    """
    cpdef Piola1Stiffness(self, double stretch, double stretch_rate, DataContainer dc, double curvature):
        cdef double P,A
        
        Me = self.material.Piola1Stiffness_1d(stretch,stretch_rate,dc.dt);
        
        P = Me['Pe']*self.area;                     # Piola stress * area
        A = Me['Ae']*self.area;                     # d(Piola stress)/d(stretch) * area
        M = Me['Ae']*self.I11/stretch*curvature     # Bending moment
        G = Me['Ae']*self.I11/stretch               # d(Bending moment)/d(curvature)
        
        return {'P':P,'A':A,'M':M,'G':G}
    """

    cpdef CauchyStiffness(self,np.ndarray[double,ndim=1] eps, np.ndarray[double,ndim=1] curvature, \
                               np.ndarray[double,ndim=2] lam, double dt,int ndim):
        cdef np.ndarray[double,ndim=1] P,M
        cdef np.ndarray[double,ndim=2] dPdF, dMdK

        # Material tangents
        dPdF = np.diag(np.array([self.material.E *self.area, \
                                 self.material.mu*self.area, \
                                 self.material.mu*self.area]))
        dPdF = lam.T.dot(dPdF.dot(lam))
        
        dMdK = np.diag(np.array([self.material.mu*(self.I11+self.I22),\
                                 self.material.E *self.I11,\
                                 self.material.E *self.I22]))
        dMdK = lam.T.dot(dMdK.dot(lam))

        # Force and moment vector
        P = dPdF.dot(eps)
        M = dMdK.dot(curvature)

        return {'P':P,'M':M,'dPdF':dPdF,'dMdK':dMdK}
