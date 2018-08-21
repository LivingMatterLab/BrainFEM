from Material cimport *
cimport numpy as np

cdef class MAT_OGDEN(Material):
	cdef public double E, nu, mu, lamb, beta, D,K
	cdef public list mu_i, alpha_i
	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] Fe, int ndim)
	cpdef Piola1Stiffness_1d(self,double stretche, double stretche_rate = *, double dt = *)
	cpdef CauchyStiffness_1d(self,double stretche, double stretche_rate = *, double dt = *)

	cpdef GetCoefficients2D(self, np.ndarray[double,ndim=1] eigVal, \
	                            np.ndarray[double,ndim=1] IC,\
	                            np.ndarray[double,ndim=1] omega, \
	                            np.ndarray[double,ndim=1] domega, \
	                            np.ndarray[double,ndim=1] ddomega, \
	                            int ndim)

	cpdef GetCoefficients3D(self, np.ndarray[double,ndim=1] eigVal, \
	                            np.ndarray[double,ndim=1] IC,\
	                            np.ndarray[double,ndim=1] omega, \
	                            np.ndarray[double,ndim=1] domega, \
	                            np.ndarray[double,ndim=1] ddomega, \
	                            int ndim)

