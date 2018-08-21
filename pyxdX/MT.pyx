# -*- coding: utf-8 -*-
cdef class MT(object):
	def __init__(self,n0,n1,n2,state=-1,timeToNextEvent=1.e9):
		self.n0=n0;self.n1=n1;self.n2=n2;
		self.e0=-1;self.e1=-1;self.e2=-1;
		self.eGC = -1;
		self.state = state
		self.timeToNextEvent = timeToNextEvent
		self.mtMinus = None; self.mtPlus = None;
		self.mechanism = None;

	def __str__(self):
		return "MT: \t [n0,n1,n2] = [" + str(self.n0) + "," + str(self.n1)+ "," + str(self.n2)+"]" +\
				  " \t [e0,e1,e2] = [" + str(self.e0) + "," + str(self.e1)+ "," + str(self.e2)+"]" +\
				  " \t eGC = " + str(self.eGC)+\
				  " \t MTstate = " + str(self.state) + "\t timeToNextEvent = " + str(self.timeToNextEvent)
