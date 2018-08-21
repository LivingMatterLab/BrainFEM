# -*- coding: utf-8 -*-
cimport numpy as np
import numpy as np

cdef class Amplitude:
	def __init__(self, t, a,dadt=None):
		self.time     	 = t
		self.amplitude   = a
		self.derivative  = dadt

	cpdef Get(self,double t):
		return np.interp(t, self.time, self.amplitude)

	cpdef GetDerivative(self,double t):
		return np.interp(t, self.time, self.derivative)