# -*- coding: utf-8 -*-\
# Bar element that generates force, assumed to be incompressible
from PBAR cimport *
from DataContainer cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PBAR1(PBAR):
    def __init__(self,material,propProp, amplitude):
        super(PBAR1,self).__init__(material,propProp);
        
        self.force  = propProp['force'];
        self.amplitude = amplitude

    def __str__(self):
        return super(PBAR1,self).__str__() + "\t Force = " + str(self.force) 

    cpdef Piola1Stiffness(self, double stretch, double stretch_rate, DataContainer dc):
        cdef double P
        
        Me = super(PBAR1,self).Piola1Stiffness(stretch,stretch_rate,dc);
        P = Me['P']-self.force*self.amplitude.Get(dc.time);         # Motor force subtracted from internal force, such that bar will lengthen
        
        return {'P':P,'A':Me['A']}

    cpdef CauchyStiffness(self, double stretch, double stretch_rate, DataContainer dc):
        cdef double sig

        Me = super(PBAR1,self).CauchyStiffness(stretch,stretch_rate,dc);
        
        sig = Me['sig']-self.force*self.amplitude.Get(dc.time); # Motor force subtracted from internal force, such that bar will lengthen
        
        return {'sig':sig,'C':Me['C']}