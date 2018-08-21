from Property cimport *
from Material cimport *
cimport numpy as np

cdef class PSOLID12(Property):
	cdef Material material
	cdef double th_rate
	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim)

