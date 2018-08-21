from DataContainer cimport *
from Element cimport *
cimport numpy as np

cdef class CBAR(Element):
	cdef public list stretch0ID, stretchID
	cdef public double restLength
	
	cpdef InitializeData(self,DataContainer dc)
	
	cpdef RotationMatrix(self)
	cpdef BuildElementMatrix(self,DataContainer dc)
	cpdef BuildInternalForceVector(self, DataContainer dc)


