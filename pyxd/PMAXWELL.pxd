from Property cimport *
from Material cimport *
cimport numpy as np

cdef class PMAXWELL(Property):
	cdef Material material0, material1
	cdef double eta
	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim)
	cpdef UpdateFg(self, np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fv0,double dt,int ndim)
    

