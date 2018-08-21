from Mechanism cimport *
from Element cimport *
from ModelContainer cimport *
from Solver cimport *

cdef class MECH_EL02(Mechanism):
	cdef public double tCrea, tDest, maxInitStretch, maxStretch
	cpdef Initialize(self, Element el, ModelContainer mc)
	cpdef Apply(self, Element el, ModelContainer mc, Solver s)
