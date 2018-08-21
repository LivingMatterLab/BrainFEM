from DataContainer cimport *
from Property cimport *
cimport numpy as np
cdef class Element:
	cdef public int localID
	cdef public str type
	cdef public list nodes, datKID, datRID
	cdef public np.ndarray direction
	cdef public Property property
	cpdef Loc(self)
	cpdef DofID(self)
	cpdef Dof0(self, DataContainer dc)
	cpdef Dof(self, DataContainer dc)

