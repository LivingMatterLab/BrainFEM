from DataContainer cimport *
from ModelContainer cimport *
cdef class Solver:
	cdef public DataContainer dc
	cdef public ModelContainer mc
	cdef public str type
	cdef public bint doParallel
	cdef public int nProc
	cpdef Solve(self)
