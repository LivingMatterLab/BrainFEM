from PBAR cimport *
from DataContainer cimport *
from Amplitude cimport *
cimport numpy as np

cdef class PBAR1(PBAR):
	cdef public double force
	cdef public Amplitude amplitude
	cpdef Piola1Stiffness(self, double stretch, double stretch_rate, DataContainer dc)
	cpdef CauchyStiffness(self, double stretch, double stretch_rate, DataContainer dc)

