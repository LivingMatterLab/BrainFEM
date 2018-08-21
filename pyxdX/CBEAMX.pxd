from CBEAM cimport *
from CBARX cimport State
from DataContainer cimport *
from Mechanism cimport *

cdef class CBEAMX(CBEAM):
	cdef public int plusID, minusID
	cdef public list dummyNodes
	cdef public double timeToNextEvent
	cdef public State state
	cdef public Mechanism mechanism

	cpdef double CurrentLength(self, DataContainer dc)