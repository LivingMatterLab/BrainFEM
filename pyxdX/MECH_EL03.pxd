from Mechanism cimport *
from Element cimport *
from ModelContainer cimport *
from Solver cimport *

cdef class MECH_EL03(Mechanism):
	cdef public double tCont, tDest, tCrea, minInitStretch, minStretch, maxInitStretch, maxStretch, activeStretch
	cpdef Initialize(self, Element el, ModelContainer mc)
	cpdef Apply(self, Element el, ModelContainer mc, Solver s)
