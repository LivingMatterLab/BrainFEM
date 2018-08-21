from Mechanism cimport *
from Element cimport *
from ModelContainer cimport *
from Solver cimport *

cdef class MECH_EL05(Mechanism):
	cdef public double tMove, minStretch, maxStretch, activeStretch
	cpdef Initialize(self, Element el, ModelContainer mc)
	cpdef Apply(self, Element el, ModelContainer mc, Solver s)
