# -*- coding: utf-8 -*-\
# Bar element, assumed to be incompressible
from Property cimport *
from DataContainer cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PBAR(Property):
    def __init__(self,material,propProp):
        super(PBAR,self).__init__('PBAR');
        
        self.area  = propProp['area'];
        self.material = material;

    def __str__(self):
        return super(PBAR,self).__str__() + "\t Area = " + str(self.area) + "\t material = " + str(self.material.localID)

    cpdef Piola1Stiffness(self, double stretch, double stretch_rate, DataContainer dc):
        cdef double P,A
        
        Me = self.material.Piola1Stiffness_1d(stretch,stretch_rate,dc.dt);
        
        P = Me['Pe']*self.area;
        A = Me['Ae']*self.area;
        
        return {'P':P,'A':A}

    cpdef CauchyStiffness(self, double stretch, double stretch_rate, DataContainer dc):
        cdef double sig,C
        
        Me = self.material.CauchyStiffness_1d(stretch,stretch_rate,dc.dt);
        
        sig = Me['sige']*self.area/stretch;
        C   = Me['Ce']*self.area/stretch;
        
        return {'sig':sig,'C':C}

    
    cpdef Piola1Stiffness0(self, double stretch, double stretch_rate, DataContainer dc):
        cdef double P,A
        
        Me = self.material.Piola1Stiffness_1d(stretch,stretch_rate,dc.dt0);
        
        P = Me['Pe']*self.area;
        A = Me['Ae']*self.area;
        
        return {'P':P,'A':A}
