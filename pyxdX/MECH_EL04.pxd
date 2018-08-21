from Mechanism cimport *
from Element cimport *
from ModelContainer cimport *
from Solver cimport *

cdef class MECH_EL04(Mechanism):
	cdef public double tCont, tDest, tCrea, minInitStretch, minStretch, maxInitStretch, maxStretch, activeStretch
	cpdef Initialize(self, Element el, ModelContainer mc)
	cpdef Apply(self, Element el, ModelContainer mc, Solver s)


	cpdef getTCrea(self, Element el, ModelContainer mc, Solver s)


cdef class MECH_EL041(MECH_EL04):
	cdef public double alpha, mtC, baseNMT

	cpdef getTCrea(self, Element el, ModelContainer mc, Solver s)
