# -*- coding: utf-8 -*-
cdef class X:
	def __init__(self,var):
		self.var = var
		print "You just created an X object with input: " , self.var