# -*- coding: utf-8 -*-
# CBARX2 is same as CBARX, but assumes that model is one dimensional
# So, it only computes element matrix and element vector along global X-direction
from CBARX cimport *
cimport numpy as np
import numpy as np
cimport ElementHelper as eh

cdef class CBARX2(CBARX):
	def __init__(self,nodes,elementProperty,state=-1,timeToNextEvent=1.e9, restLength=0.0):
		super(CBARX2,self).__init__(nodes,elementProperty,state,timeToNextEvent, restLength);

	cpdef BuildElementMatrix(self,DataContainer dc):
		cdef double l, stretch, stretch_rate, nn,Inn
		cdef np.ndarray[double,ndim=1] ex0, ex
		cdef np.ndarray[double,ndim=2] mat,T,mat2
		cdef dict M

		# If inactive MT, don't compute
		if self.state == MicrotubuleInactive:
			mat = 100000.*np.array([[1., -1.],[-1., 1.]])
			return np.asarray(mat)


		# Compute stretch as ratio of direction vectors
		if (self.state == Microtubule):
			l = self.nodes[1].x-self.nodes[0].x + \
				self.nodes[1].Dof(dc)[0]-self.nodes[0].Dof(dc)[0]
		else:
			ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
			ex  = ex0
			ex[0] += self.nodes[1].Dof(dc)[0]-self.nodes[0].Dof(dc)[0]
			l   = np.linalg.norm(ex);
		stretch = l/self.restLength;
		stretch_rate = (stretch-eh.getStretch0(self,0,dc))/dc.dt;

		if stretch<1e-6:
			print self
			print self.nodes[0]
			print self.nodes[0].Dof(dc)
			print self.nodes[1]
			print self.nodes[1].Dof(dc)

		# Compute normal force and tangent
		M = self.property.Piola1Stiffness(stretch,stretch_rate,dc);

		if (self.state == Microtubule):
			mat = M['A']/self.restLength * np.array([[1., -1.],[-1., 1.]])
		else:
			# Compute element matrix
			nn = ex[0]*ex[0]/(l*l)
			Inn = 1-nn

			mat = M['A']/self.restLength * np.array([[ nn, -1* nn],[-1* nn, nn]]) + \
			      M['P']/l  * np.array([[Inn, -1*Inn],[-1*Inn,Inn]])
		
		return np.asarray(mat)
		
	cpdef BuildInternalForceVector(self,DataContainer dc):
		cdef double l, stretch, stretch_rate
		cdef np.ndarray[double,ndim=1] ex0, ex, vec
		cdef dict M

		# Compute stretch as ratio of direction vectors
		ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
		ex  = ex0
		ex[0] += self.nodes[1].Dof(dc)[0]-self.nodes[0].Dof(dc)[0]
		l   = np.linalg.norm(ex);
		stretch = l/self.restLength;
		eh.setStretch(self,0,stretch,dc)
		stretch_rate = (stretch-eh.getStretch0(self,0,dc))/dc.dt;

		# Compute normal force and tangent
		M = self.property.Piola1Stiffness(stretch,stretch_rate,dc);

		# Compute residual force vector
		vec = M['P']/l* np.array([-ex[0],ex[0]])

		return vec