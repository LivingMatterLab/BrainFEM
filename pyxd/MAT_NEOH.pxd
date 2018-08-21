from Material cimport *
cimport numpy as np

cdef class MAT_NEOH(Material):
	cdef public double E, nu, mu, lamb, D 
	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] Fe, int ndim)
	cpdef Piola1Stiffness_1d(self,double stretche, double stretche_rate = *, double dt = *)
	cpdef CauchyStiffness_1d(self,double stretche, double stretche_rate = *, double dt = *)

