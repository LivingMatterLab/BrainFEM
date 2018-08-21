from Property cimport *
from Material cimport *
from DataContainer cimport *
cimport numpy as np

cdef class PBAR(Property):
	cdef public Material material
	cdef public double area
	cpdef Piola1Stiffness(self, double stretch, double stretch_rate, DataContainer dc)
	cpdef Piola1Stiffness0(self, double stretch, double stretch_rate, DataContainer dc)
	cpdef CauchyStiffness(self, double stretch, double stretch_rate, DataContainer dc)

