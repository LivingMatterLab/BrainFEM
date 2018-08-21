from CBAR cimport *
from DataContainer cimport *
from Mechanism cimport *

cdef enum State:
	NoState       = -1
	Actin         =  15
	Microtubule   =  10 
	MicrotubuleInactive   =  -10 
	MotorAttachedSlack =  1
	MotorAttachedTaut  = 2
	MotorReleased =  3
	MotorJustAttachedSlack = 4
	MotorJustAttachedTaut  = 5
	MotorJustReleased = 6

cdef class CBARX(CBAR):
	cdef public int plusID, minusID
	cdef public list dummyNodes
	cdef public double timeToNextEvent
	cdef public State state
	cdef public Mechanism mechanism

	cpdef double CurrentLength(self, DataContainer dc)