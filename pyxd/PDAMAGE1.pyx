# -*- coding: utf-8 -*-\
# Anisotropic prescribed growth: Fg = I + (theta-1)*e1*e1
from Property cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PDAMAGE1(Property):
    def __init__(self,material,propProp):
        super(PDAMAGE1,self).__init__('PDAMAGE1');
        
        self.aL  = propProp['aL'];
        self.bL  = propProp['bL'];
        self.aA  = propProp['aA'];
        self.bA  = propProp['bA'];

        self.material = material;

    def __str__(self):
        return super(PDAMAGE1,self).__str__() + "\t aL = " + str(self.aL) \
                                              + "\t bL = " + str(self.bL) \
                                              + "\t aA = " + str(self.aA) \
                                              + "\t bA = " + str(self.bA) \
                                              + "\t material = " + str(self.material.localID)

    cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim):
        cdef int i,j,k,l
        cdef double stretch
        cdef np.ndarray[double, ndim=1] e1
        cdef np.ndarray[double, ndim=2] Fg, Fginv, P, Pe
        cdef np.ndarray[double, ndim=4] A, Ae

        stretch = F[0,0]
        e1 = np.zeros(ndim); e1[0] = 1.;
        Fg = Fg0;

        # Initialize damage as 0,0 component of Fg
        if Fg[0,0]==1.:
            Fg[0,0]=0;

        # Compute damage based on stretch
        d = 1./(1.+np.exp(-self.aA*(stretch-self.aL)))

        # Check if damage is larger than current damage.
        # If so, update damage. Otherwise damage remains same (cannot decrease)
        if d>Fg[0,0]:
            Fg[0,0] = d;
            dd      = self.aA*(d-d**2);
        else:
            d = Fg[0,0]
            dd = 0;

        Me = self.material.Piola1Stiffness(F,ndim);
        

        Pe = Me['Pe']
        P  = (1-d)*Pe


        Ae = Me['Ae']
        A  = (1-d)*Ae -dd/stretch*np.einsum('ij,kl->ijkl',Pe,np.outer(F.dot(e1),e1))
    
        return {'P':P,'A':A,'Fg':Fg}
