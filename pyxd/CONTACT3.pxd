from DataContainer cimport *
cimport numpy as np

cdef class CONTACT:
	cdef public double beta1
	cdef list slaveNodes, masterNodes
	cdef public int numSlaveNodes, numMasterNodes
	cdef public double penaltyParameter
	cdef public list datRID, datKID

	cpdef BuildMatrix(self, DataContainer dc)
	cpdef BuildVector(self, DataContainer dc)

	
	cpdef checkSlaveContact(self, np.ndarray[double,ndim=1] sLoc, list mLoc, list lenElem, list tangVec, list normVec)

	cpdef getHermiteG(self, double g,double alpha,\
		           np.ndarray[double,ndim=1] x1,\
		           np.ndarray[double,ndim=1] x2,\
		           np.ndarray[double,ndim=1] x3,\
		           np.ndarray[double,ndim=1] x4)

	cpdef getHermiteDepth(self, double g,double alpha,\
		           np.ndarray[double,ndim=1] x1,\
		           np.ndarray[double,ndim=1] x2,\
		           np.ndarray[double,ndim=1] x3,\
		           np.ndarray[double,ndim=1] x4)

	cpdef getPenalty(self,double penDepth)
	
