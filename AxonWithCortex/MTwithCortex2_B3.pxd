from MTwithCortex_B3 cimport *

cdef class MTwithCortex2_B3(MTwithCortex_B3):
	cpdef BuildModel(self,object p)

	cpdef CreateMT(self, object p, dict geom, int nCount, int eCount)
	cpdef CreateMTCrosslinks(self, object p, dict geom, int nCount, int eCount)

	cpdef CreateActin(self, object p, dict geom, int nCount, int eCount);
	cpdef CreateActinCrosslinks(self,object p,dict geom, list nidRing0, list nidRing1, int nCount, int eCount)
	cpdef CreateAcMTCrosslinks(self, object p, dict geomMT, dict geomAc,\
	                           dict geomCL_Ac_Mt, int nCount, int eCount)
