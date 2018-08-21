from DataContainer cimport *
from Element cimport *
cimport numpy as np

cdef class CBEAM(Element):
	cdef public list stretch0ID, stretchID, curvature0ID, curvatureID, theta0ID, thetaID
	cdef int nip
	cdef public double restLength
	cdef np.ndarray T
	
	cpdef InitializeData(self,DataContainer dc)
	
	cpdef RotationMatrix(self)
	cpdef BuildElementMatrix(self,DataContainer dc)
	cpdef BuildInternalForceVector(self, DataContainer dc)

	cpdef T_thw(self,np.ndarray[double,ndim=1] th_vec,int ndim)
	cpdef T_wth(self,np.ndarray[double,ndim=1] th_vec,int ndim)
	cpdef UpdateTheta(self,int ip,np.ndarray[double,ndim=1] dw, DataContainer dc,int ndim)


