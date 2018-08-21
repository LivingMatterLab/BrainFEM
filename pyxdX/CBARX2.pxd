from CBARX cimport *
from DataContainer cimport *

cdef class CBARX2(CBARX):
	cpdef BuildElementMatrix(self,DataContainer dc)
	cpdef BuildInternalForceVector(self, DataContainer dc)