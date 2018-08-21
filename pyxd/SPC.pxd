from Node cimport *
from Amplitude cimport *

cdef class SPC:
	cdef public Node node
	cdef public list dofIDs
	cdef public double value
	cdef public Amplitude amplitude

	cpdef Get(self, double time, double dt)

