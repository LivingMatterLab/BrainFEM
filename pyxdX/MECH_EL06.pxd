from Mechanism cimport *
from Element cimport *
from ModelContainer cimport *
from Solver cimport *

cdef class MECH_EL06(Mechanism):
	cdef public double tDest0, f_beta, maxInitStretch, maxStretch
	cdef double tCountDown
	cpdef Initialize(self, Element el, ModelContainer mc)
	cpdef Apply(self, Element el, ModelContainer mc, Solver s)
	
	cpdef CdfStatBreakage(self, double time, double force)

	cpdef PdfDynBreakage(self, double force, double force_rate)
	cpdef CdfDynBreakage(self, double force, double force_rate)
