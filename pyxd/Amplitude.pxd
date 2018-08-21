cimport numpy as np
cdef class Amplitude:
	cdef public np.ndarray time, amplitude, derivative;
	cpdef Get(self, double t)
	cpdef GetDerivative(self, double t)