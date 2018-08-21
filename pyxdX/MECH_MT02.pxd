from MECH_MT01 cimport *
from MT cimport *
from ModelContainer cimport *
from Solver cimport *

cdef class MECH_MT02(MECH_MT01):
	cdef public double thresDist
	cdef public bint stabDistalMT
	cpdef Initialize(self, MT mt, ModelContainer mc)
	cpdef Apply(self, MT mt, ModelContainer mc, Solver s)

	cpdef InitializeGCNode(self, MT mt, ModelContainer mc)
	cpdef ApplyGCNode(self, MT mt, ModelContainer mc, Solver s)
