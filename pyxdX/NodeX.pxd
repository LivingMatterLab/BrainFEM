from Node cimport *
from Element cimport *
from MPC cimport *

cdef class NodeX(Node):
	cdef public Element elPlus, elMinus
	cdef public MPC mpc
	cdef public double curvature0

	cpdef InitializeCurvature(self,DataContainer dc)
