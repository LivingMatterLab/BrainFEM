# -*- coding: utf-8 -*-
cdef class SPC:
	def __init__(self, node, dofIDs,value, amplitude):
		self.node = node
		self.dofIDs = dofIDs
		self.value = value
		self.amplitude = amplitude


	def __str__(self):
		return "SPC: \t Node " + str(self.node.localID)+"  \t dofIDs = [" + \
				', '.join(str(x) for x in self.dofIDs)+"]" +  \
				"\t value = "+ str(self.value)

	cpdef Get(self, double time, double dt):
		return self.value * (self.amplitude.Get(time) - self.amplitude.Get(time-dt))
