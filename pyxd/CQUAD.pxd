from DataContainer cimport *
from Element cimport *
cimport numpy as np

cdef class CQUAD(Element):
	cdef np.ndarray T
	cdef public list Fg0ID, FgID
	cdef bint reducedIntegration


	cpdef InitializeData(self,DataContainer dc)
	
	cpdef RotationMatrix(self)
	cpdef BuildElementMatrix(self,DataContainer dc)
	cpdef BuildInternalForceVector(self, DataContainer dc)


