from Solver cimport *
from ModelContainer cimport *
cdef class SOL2(Solver):
	cdef public int maxIter, maxIterInc, maxStep
	cdef public double tEnd, dt0, dtMin, dtMax
	cdef public double[:] datK, datR
	cdef public bint plotOutput
	cpdef Solve(self)

	cpdef StiffnessData(self, list eids, list pids)
	cpdef StiffnessDataCont(self,list cids)
	cpdef ResidualDataExt(self, list lids, list pids, double time)
	cpdef ResidualDataCont(self, list cids, int it)
	cpdef ResidualDataInt(self, list eids)
	cpdef MPCMatrix(self, ModelContainer mc)
	cpdef RowColumnVectors(self)


