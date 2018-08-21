# -*- coding: utf-8 -*-
cdef class LOAD:
	def __init__(self, node, loadVec, amplitude):
		self.node = node
		self.loadVec = loadVec
		self.amplitude = amplitude

	def __str__(self):
		return "LOAD: \t Node " + str(self.node.localID)+"  \t loadVec = ["+ ', '.join(str(x) for x in self.loadVec)+"]"

	cpdef Get(self, double time):
		return self.loadVec*self.amplitude.Get(time)
