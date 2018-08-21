cimport numpy as np
from Element cimport *
from DataContainer cimport *
from Solver cimport *

cpdef LocalDeformation_B3(double r, np.ndarray[double,ndim=1] r0, np.ndarray[double,ndim=1] v)

cpdef LocalDeformation_Q2(double r, double s, np.ndarray[double,ndim=1] r0, np.ndarray[double,ndim=1] v)

cpdef LocalDeformation_H3(double r, double s, double t, np.ndarray[double,ndim=1] r0, np.ndarray[double,ndim=1] v)

cpdef Jacobian_B3(double r,  np.ndarray[double,ndim=1] r0)

cpdef Jacobian_Q2(double r, double s, np.ndarray[double,ndim=1] r0)

cpdef Jacobian_H3(double r, double s, double t, np.ndarray[double,ndim=1] r0)

cpdef GaussPoints(int nPoints)

cpdef SkewToMat(np.ndarray[double,ndim=1] vec)

cpdef getFg0(Element el, int ip, DataContainer dc, int ndim)
cpdef getFg(Element el, int ip, DataContainer dc, int ndim)

cpdef getCurvature0(Element el, int ip, DataContainer dc, int ndim)
cpdef getCurvature(Element el, int ip, DataContainer dc, int ndim)

cpdef getTheta0(Element el, int ip, DataContainer dc, int ndim)
cpdef getTheta(Element el, int ip, DataContainer dc, int ndim)

cpdef getStretch0(Element el, int ip, DataContainer dc)
cpdef getStretch(Element el, int ip, DataContainer dc)

cpdef getNc0(Element el, int ip, DataContainer dc)
cpdef getNc(Element el, int ip, DataContainer dc)

cpdef getJg(Element el, int ip, DataContainer dc, int ndim)

cpdef setFg(Element el, int ip, np.ndarray[double,ndim=2] Fg, DataContainer dc, int ndim)
cpdef setCurvature(Element el, int ip, np.ndarray[double,ndim=1] curvature, DataContainer dc, int ndim)
cpdef setTheta(Element el, int ip, np.ndarray[double,ndim=1] theta, DataContainer dc, int ndim)
cpdef setStretch(Element el, int ip, double stretch, DataContainer dc)
cpdef setNc(Element el, int ip, double nc, DataContainer dc)

cpdef UpdateConnectivities(Element el, Solver s)