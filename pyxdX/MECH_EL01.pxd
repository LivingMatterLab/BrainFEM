from Mechanism cimport *
from Element cimport *
from ModelContainer cimport *
from Solver cimport *

cdef class MECH_EL01(Mechanism):
	cdef public double tCrea, tDest
	cpdef Apply(self, Element el, ModelContainer mc, Solver s)
