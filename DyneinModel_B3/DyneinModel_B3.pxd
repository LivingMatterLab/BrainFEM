from DataContainer cimport *
from ModelContainer cimport *
from Solver cimport *
from Node cimport *

cdef class DyneinModel_B3(ModelContainer):
	cdef public list storageNodes, MTs, elRestore, mtRestore
	cdef Node loadNode
	cdef int eidCL0, eidCL1
	cdef public int optionDynein
	cdef public double lAxon
	cpdef BuildModel(self,object p)
	cpdef getGeometry(self,object p)
	cpdef WriteStepOutput(self, Solver s)
	cpdef UpdateModel(self, Solver s)
	cpdef RestoreModel(self, Solver s)
