from Property cimport *
from Material cimport *
from DataContainer cimport *
from Amplitude cimport *
cimport numpy as np

cdef class PSOLID71(Property):
	cdef public Material material
	cdef public double Grho, advSpeed, kth1, kth2, alpha
	cdef public Amplitude amplRhoSource, amplAdvThreshold

	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double Rho, double nc0, double dt, int ndim)
	cpdef Diffusivity(self,np.ndarray[double,ndim=2] F, double Rho, np.ndarray[double,ndim=1] gradRho, int ndim)
	cpdef SourceRho(self,np.ndarray[double,ndim=2] F, double Rho, np.ndarray[double,ndim=1] gradRho, DataContainer dc, int ndim)
	cpdef UpdateNc(self,double nc0, double dt)
	cpdef UpdateTh1(self,double nc, double Rho)
	cpdef UpdateTh2(self,double nc, double Rho)

