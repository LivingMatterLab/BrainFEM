from ModelContainer cimport *

cdef class DataContainer:
	cdef public double[:] dof, dof0
	cdef public double[:] Rp, Rp0
	cdef public double time, dt, dt0
	cdef public int step
	cdef public list Fg0, Fg, curvature0, curvature, theta0, theta,stretch0, stretch, nc0, nc, OutputAnimation
	cdef public double[:] Fg0_, Fg_, curvature0_, curvature_,theta0_, theta_, stretch0_, stretch_, nc0_, nc_
	cdef public list countFg0, countFg, countCurvature0, countCurvature, countTheta0, countTheta, countStretch0, countStretch, countNc0, countNc

	cpdef InitializeData(self, ModelContainer mc)
	cpdef ShareData(self)
	cpdef AddToData(self, listToAddTo, countToAddTo, dataToAdd)

