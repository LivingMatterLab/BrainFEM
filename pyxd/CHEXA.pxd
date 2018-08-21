from DataContainer cimport *
from Element cimport *
cimport numpy as np

cdef class CHEXA(Element):
	cdef np.ndarray T
	cdef public list Fg0ID, FgID

	cpdef InitializeData(self,DataContainer dc)
	cpdef RotationMatrix(self)
	cpdef BuildElementMatrix(self,DataContainer dc)
	cpdef BuildInternalForceVector(self, DataContainer dc)


