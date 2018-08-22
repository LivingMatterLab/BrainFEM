from ModelContainer cimport *
from Solver cimport *
cdef class BrainDevCirc_Q2(ModelContainer):
	cpdef BuildModel(self,object p)
	
	cpdef GetElemDir(self, double x, double y, str option)
	cpdef ReadAbaqusFile(self,str filename)

	cpdef WriteStepOutput(self, Solver s)
