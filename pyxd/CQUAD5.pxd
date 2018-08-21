from DataContainer cimport *
from CQUAD cimport *
cimport numpy as np

cdef class CQUAD5(CQUAD):
	cdef np.ndarray R0

	cpdef InitializeData(self,DataContainer dc)

	cpdef ComputeR0(self)
	cpdef BuildElementMatrix(self,DataContainer dc)
	cpdef BuildInternalForceVector(self, DataContainer dc)


