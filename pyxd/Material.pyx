# -*- coding: utf-8 -*-
cdef class Material(object):
	def __init__(self,materialType):
		self.localID = -1
		self.type = materialType

	def __str__(self):
		return "Material " + str(self.localID) +":  \t type = " + self.type