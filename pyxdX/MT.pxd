from Mechanism cimport *

cdef enum MTstate:
	NoState = -1
	Polymerizing = 0
	Stationary =  1 
	Depolymerizing =  2 

cdef class MT:
	cdef public int n0,n1,n2,e0,e1,e2,eGC
	cdef public MTstate state
	cdef public double timeToNextEvent
	cdef public MT mtMinus, mtPlus
	cdef public Mechanism mechanism