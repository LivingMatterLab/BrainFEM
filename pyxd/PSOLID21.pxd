from Property cimport *
from Material cimport *
cimport numpy as np

cdef class PSOLID21(Property):
	cdef Material material
	cdef double Gs, stretchCrit
	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim)
	cpdef UpdateFg(self, np.ndarray[double,ndim=2] F,double th0,double dt)

