# -*- coding: utf-8 -*-
cdef class Solver(object):
	def __init__(self, solType, mc, dc,doParallel,nProc=4):
		self.type = solType
		self.mc = mc
		self.dc = dc
		self.doParallel = doParallel
		self.nProc = nProc


	cpdef Solve(self):
		print "ERROR: BuildModel function not implemented in subclass"


