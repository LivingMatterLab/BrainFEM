from DataContainer cimport *
from Element cimport *
cimport numpy as np

cdef class CQUAD7(Element):
	cdef np.ndarray T
	cdef public list Fg0ID, FgID, nc0ID, ncID


	cpdef InitializeData(self,DataContainer dc)
	
	cpdef RotationMatrix(self)
	cpdef BuildElementMatrix(self,DataContainer dc)
	cpdef BuildInternalForceVector(self, DataContainer dc)

	cpdef getVolume(self, DataContainer dc)
	cpdef getDensity(self, DataContainer dc)


