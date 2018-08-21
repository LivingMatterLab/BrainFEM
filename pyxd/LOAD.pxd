from Node cimport *
from Amplitude cimport *

cdef class LOAD:
	cdef public Node node
	cdef public np.ndarray loadVec
	cdef public list datRID
	cdef public Amplitude amplitude

	cpdef Get(self, double time)


