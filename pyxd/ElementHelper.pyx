# -*- coding: utf-8 -*-
# ============================================ #
# This file contains general FUNCTIONS that 
# return variables that are very standard in 
# FEM, and hence are required by most elements
# that will be implemented. These are e.g. 
# shape functions, their derivatives, jacobians,
# etc.
# -------------------------------------------- #
cimport numpy as np
import numpy as np
from Element cimport *
from DataContainer cimport *

cpdef LocalDeformation_B3(double r,np.ndarray[double,ndim=1] r0,np.ndarray[double,ndim=1] v):
	cdef np.ndarray[double, ndim=1] r_c, dNx,F,DL
	cdef dict jac

	# First compute Jacobian, its inverse, and Ndr
	jac = Jacobian_B3(r,r0)

	# Derivatives of N wrt global coordinates
	dNx = jac['dNr']*jac['Jinv'];

	# Deformation gradient
	r_c = r0+v;
	F = np.reshape(r_c,(2,3)).T.dot(dNx)

	return {'F':F,'J':jac['J'],'dNx':dNx,'N':jac['N']}

cpdef LocalDeformation_Q2(double r,double s,np.ndarray[double,ndim=1] r0,np.ndarray[double,ndim=1] v):
	cdef np.ndarray[double, ndim=1] r_c
	cdef np.ndarray[double, ndim=2] dNx, F
	cdef dict jac

	# First compute Jacobian, its inverse, and Ndr
	jac = Jacobian_Q2(r,s,r0)

	# Derivatives of N wrt global coordinates
	dNx = jac['dNr'].dot(jac['Jinv']);

	# Deformation gradient
	r_c = r0+v;
	F = np.reshape(r_c,(4,2)).T.dot(dNx)

	return {'F':F,'J':jac['J'],'dNx':dNx,'N':jac['N']}

cpdef LocalDeformation_H3(double r,double s,double t,np.ndarray[double,ndim=1] r0,np.ndarray[double,ndim=1] v):
	cdef np.ndarray[double, ndim=1] r_c
	cdef np.ndarray[double, ndim=2] dNx, F
	cdef dict jac

	# First compute Jacobian, its inverse, and Ndr
	jac = Jacobian_H3(r,s,t,r0)

	# Derivatives of N wrt global coordinates
	dNx = jac['dNr'].dot(jac['Jinv']);

	# Deformation gradient
	r_c = r0+v;
	F = np.reshape(r_c,(8,3)).T.dot(dNx)

	return {'F':F,'J':jac['J'],'dNx':dNx,'N':jac['N']}

cpdef Jacobian_B3(double r, np.ndarray[double,ndim=1] r0):
	cdef np.ndarray[double, ndim=1] dNdr, dNds, N, dNr
	cdef np.ndarray[double, ndim=2] gp
	cdef double J, Jinv
	cdef int i

	gp = np.array([[-1.],[1.]])

	N = (1+gp[:,0]*r)/2;
	dNr = np.empty(2)
	for i in range(2):
		dNr[i] = ( gp[i,0])/2.0

	J = (r0[3]-r0[0])/2.
	Jinv = 1/J;

	return {'J':J,'Jinv':Jinv,'dNr':dNr,'N':N}

cpdef Jacobian_Q2(double r,double s, np.ndarray[double,ndim=1] r0):
	cdef np.ndarray[double, ndim=1] dNdr, dNds, N
	cdef np.ndarray[double, ndim=2] gp, dNr, J, Jinv
	cdef int i

	gp = np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])

	N = (1+gp[:,0]*r)*(1+gp[:,1]*s)/4;

	dNdr = np.empty(4)
	dNds = np.empty(4)
	for i in range(4):
		dNdr[i] = ( gp[i,0] 		* (1+gp[i,1]*s)	)/4.0
		dNds[i] = ( (1+gp[i,0]*r) 	* gp[i,1]		)/4.0

	dNr = np.squeeze([dNdr,dNds]).T
	J = np.reshape(r0,(4,2)).T.dot(dNr)
	Jinv = np.linalg.inv(J);

	return {'J':J,'Jinv':Jinv,'dNr':dNr,'N':N}

cpdef Jacobian_H3(double r,double s,double t,np.ndarray[double,ndim=1] r0):
	cdef np.ndarray[double, ndim=1] dNdr, dNds, dNdt
	cdef np.ndarray[double, ndim=2] gp, dNr, J, Jinv
	cdef int i

	gp = np.array([[-1., -1., -1.],[1., -1., -1.],[1., 1., -1.], [-1., 1., -1.],\
	               [-1., -1., 1.], [1., -1., 1.],[1., 1., 1.],[-1., 1., 1.]])

	N = (1+gp[:,0]*r)*(1+gp[:,1]*s)*(1+gp[:,2]*t)/8.;

	dNdr = np.empty(8)
	dNds = np.empty(8)
	dNdt = np.empty(8)
	for i in range(8):
		dNdr[i] = ( gp[i,0] 		* (1+gp[i,1]*s) 	* (1+gp[i,2]*t) )/8.0
		dNds[i] = ( (1+gp[i,0]*r) 	* gp[i,1]		 	* (1+gp[i,2]*t) )/8.0
		dNdt[i] = ( (1+gp[i,0]*r)	* (1+gp[i,1]*s) 	* gp[i,2]		)/8.0

	dNr = np.squeeze([dNdr,dNds,dNdt]).T
	J = np.reshape(r0,(8,3)).T.dot(dNr)
	Jinv = np.linalg.inv(J);

	return {'J':J,'Jinv':Jinv,'dNr':dNr,'N':N}

cpdef GaussPoints(int nPoints):
	cdef np.ndarray[double, ndim=1] x, w
	if(nPoints==1):
		x = np.array([0.0]);
		w = np.array([2.0]);
	elif(nPoints==2):
		x = np.array([-0.5773502691896257,0.5773502691896257]);
		w = np.array([1.0 ,1.0]);
	elif(nPoints==3):
		x = np.array([-0.7745966692414834,0,0.7745966692414834]);
		w = np.array([0.5555555555555556,0.8888888888888888,0.5555555555555556]);
	elif(nPoints==4):
		x = np.array([-0.8611363115940526,-0.3399810435848563,0.3399810435848563,0.8611363115940526])
		w = np.array([0.3478548451374538,0.6521451548625461,0.6521451548625461,0.3478548451374538])
	else:
		x = 1/0;
		w = 1/0;
	return (x,w)

cpdef SkewToMat(np.ndarray[double,ndim=1] vec):
	return np.array([[ 0.    , vec[2],-vec[1]],\
			         [-vec[2], 0.    , vec[0]],\
			         [ vec[1], -vec[0], 0.]])

cpdef getFg0(Element el, int ip, DataContainer dc, int ndim):
	return np.array([dc.Fg0_[el.Fg0ID[ip]+i] for i in range(ndim*ndim)]).reshape(ndim,ndim)

cpdef getFg(Element el, int ip, DataContainer dc, int ndim):
	return np.array([dc.Fg_[el.FgID[ip]+i] for i in range(ndim*ndim)]).reshape(ndim,ndim)

cpdef getCurvature0(Element el, int ip, DataContainer dc, int ndim):
	return np.array([dc.curvature0_[el.curvature0ID[ip]+i] for i in range(ndim)])

cpdef getCurvature(Element el, int ip, DataContainer dc, int ndim):
	return np.array([dc.curvature_[el.curvatureID[ip]+i] for i in range(ndim)])

cpdef getTheta0(Element el, int ip, DataContainer dc, int ndim):
	return np.array([dc.theta0_[el.theta0ID[ip]+i] for i in range(ndim)])

cpdef getTheta(Element el, int ip, DataContainer dc, int ndim):
	return np.array([dc.theta_[el.thetaID[ip]+i] for i in range(ndim)])

cpdef getStretch0(Element el, int ip, DataContainer dc):
	return dc.stretch0_[el.stretch0ID[ip]]

cpdef getStretch(Element el, int ip, DataContainer dc):
	return dc.stretch_[el.stretchID[ip]]

cpdef getNc0(Element el, int ip, DataContainer dc):
	return dc.nc0_[el.nc0ID[ip]]

cpdef getNc(Element el, int ip, DataContainer dc):
	return dc.nc_[el.ncID[ip]]

cpdef getJg(Element el, int ip, DataContainer dc, int ndim):
	return np.linalg.det(getFg(el,ip,dc,ndim))

cpdef setFg(Element el, int ip, np.ndarray[double,ndim=2] Fg, DataContainer dc, int ndim):
	cdef int i,j
	for i in range(ndim):
		for j in range(ndim):
			dc.Fg_[el.FgID[ip]+ndim*i+j] = Fg[i,j]

cpdef setCurvature(Element el, int ip, np.ndarray[double,ndim=1] curvature, DataContainer dc, int ndim):
	cdef int i
	for i in range(ndim):
		dc.curvature_[el.curvatureID[ip]+i] = curvature[i]

cpdef setTheta(Element el, int ip, np.ndarray[double,ndim=1] theta, DataContainer dc, int ndim):
	cdef int i
	for i in range(ndim):
		dc.theta_[el.thetaID[ip]+i] = theta[i]

cpdef setStretch(Element el, int ip, double stretch, DataContainer dc):
	dc.stretch_[el.stretchID[ip]] = stretch

cpdef setNc(Element el, int ip, double nc, DataContainer dc):
	dc.nc_[el.ncID[ip]] = nc;

cpdef UpdateConnectivities(Element el, Solver s):
	cdef np.ndarray[np.int_t, ndim=1] dof
	cdef int i, ndof, countK
	
	# Update rowR, rowK, colK in the solver
	dof = el.DofID();
	ndof = len(dof);

	countK = el.datKID[0]
	for i in range(ndof):
		s.rowR[el.datRID[0]+i]=dof[i]
		for j in range(ndof):
			s.rowK[countK] = dof[j]
			s.colK[countK] = dof[i]
			countK +=1
