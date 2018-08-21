from DataContainer cimport *
from Element cimport *
cimport numpy as np

cdef class CBAR1(Element):
	cdef np.ndarray T
	cdef public list stretch0ID, stretchID
	
	cpdef InitializeData(self,DataContainer dc)
	
	cpdef RotationMatrix(self)
	cpdef BuildElementMatrix(self,DataContainer dc)
	cpdef BuildInternalForceVector(self, DataContainer dc)


