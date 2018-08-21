from Solver cimport *
from ModelContainer cimport *
cdef class SOL2X(Solver):
	cdef public int maxIter, maxIterInc, maxStep
	cdef public double tEnd, dt0, dtMin, dtMax, tolerance
	cdef public double[:] datK, datR
	cdef public list rowK, rowR, colK, colR
	cdef public bint plotOutput
	cpdef Solve(self)

	cpdef StiffnessData(self, list eids)
	cpdef ResidualDataExt(self, list lids, double time)
	cpdef ResidualDataInt(self, list eids)
	cpdef MPCMatrix(self, ModelContainer mc)
	cpdef RowColumnVectors(self)