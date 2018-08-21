from DataContainer cimport *

cdef class CONTACT:
	cdef public double beta1
	cdef list slaveNodes, masterNodes
	cdef public int numSlaveNodes, numMasterNodes
	cdef public double penaltyParameter
	cdef public list datRID, datKID

	cpdef BuildMatrix(self, DataContainer dc)
	cpdef BuildVector(self, DataContainer dc)

