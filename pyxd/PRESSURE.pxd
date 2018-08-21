from DataContainer cimport *
from Amplitude cimport *
cimport numpy as np

cdef class PRESSURE:
	cdef public double beta1
	cdef list nodes
	cdef public int numNodes
	cdef public double pressure
	cdef public list datRID, datKID
	cdef public Amplitude amplitude

	cpdef BuildMatrix(self, DataContainer dc)
	cpdef BuildVector(self, DataContainer dc)
	
