from DataContainer cimport *
from ModelContainer cimport *
cdef class BarTest_B3(ModelContainer):
	cpdef BuildModel(self, object p)
	cpdef TestElementMatrix(self,DataContainer dc)
	cpdef PostProcess(self, DataContainer dc)
