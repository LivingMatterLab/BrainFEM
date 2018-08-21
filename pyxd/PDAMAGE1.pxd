from Property cimport *
from Material cimport *
cimport numpy as np

cdef class PDAMAGE1(Property):
	cdef Material material
	cdef double aL,bL,aA,bA
	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim)

