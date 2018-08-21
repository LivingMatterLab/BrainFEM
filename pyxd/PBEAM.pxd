from Property cimport *
from Material cimport *
from DataContainer cimport *
cimport numpy as np

cdef class PBEAM(Property):
	cdef public Material material
	cdef public double area, I11, I22, I12
	
	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=1] F, np.ndarray[double,ndim=1] curvature,double dt,int ndim)
	

	cpdef CauchyStiffness(self,np.ndarray[double,ndim=1] eps, np.ndarray[double,ndim=1] curvature ,np.ndarray[double,ndim=2] lam, double dt,int ndim)

