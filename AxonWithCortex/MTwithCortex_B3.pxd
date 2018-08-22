from DataContainer cimport *
from ModelContainer cimport *
from Solver cimport *
from Node cimport *

cdef class MTwithCortex_B3(ModelContainer):
	cdef public list storageNodes, MTs, Acs, elRestore, mtRestore
	cdef public Node loadNodeMT, growthConeNode
	cdef int eidMTCL0, eidMTCL1, eidAcCL0, eidAcCL1, eidMTAcCL0, eidMTAcCL1
	cdef public int optionDynein
	cdef public double lAxon, currLMT
	cpdef BuildModel(self,object p)
	
	cpdef getGeometry(self,object p)
	cpdef getSimpleCortexGeometry(self,object p)
	cpdef getCortexGeometry(self,object p)
	cpdef getCL_MT_Ac(self,object p, dict geomMT, dict geomAc)

	cpdef CreateMT(self, object p, dict geom, int nCount, int eCount)
	cpdef CreateMTCrosslinks(self, object p, dict geom, int nCount, int eCount)


	cpdef CreateActin(self, object p, dict geom, int nCount, int eCount);
	cpdef CreateActinCrosslinks(self,object p,dict geom, list nidRing0, list nidRing1, int nCount, int eCount)

	cpdef CalcLMT(self, Solver s)

	cpdef WriteStepOutput(self, Solver s)
	cpdef UpdateModel(self, Solver s)
	cpdef RestoreModel(self, Solver s)
