from Mechanism cimport *
from MT cimport *
from ModelContainer cimport *
from Solver cimport *

cdef class MECH_MT01(Mechanism):
	cdef public double tMTpoly,tMTstat,tMTdepoly,MTpolyRate,MTdepolyRate
	cdef public double fracLMT, goalNMT
	cpdef Initialize(self, MT mt, ModelContainer mc)
	cpdef Apply(self, MT mt, ModelContainer mc, Solver s)

	cpdef ChangeMTState(self, MT mt, ModelContainer mc, Solver s)
	cpdef GetNumMTinCS(self, MT mt, ModelContainer mc, Solver s)
