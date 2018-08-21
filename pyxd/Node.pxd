from DataContainer cimport *
cimport numpy as np

cdef class Node:
	cdef public int localID
	cdef public double x,y,z
	cdef public list loc, dofID
	cpdef Dof(self, DataContainer dc)
	cpdef Dof0(self, DataContainer dc)
	cpdef DispDof(self, DataContainer dc)
	cpdef DispDof0(self, DataContainer dc)
