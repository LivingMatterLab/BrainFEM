# -*- coding: utf-8 -*-
# Bar element for total lagrangian method
from Element cimport *
from math import *
cimport ElementHelper as eh

cimport numpy as np
import numpy as np
import scipy.linalg as spl

cdef class CBAR1(Element):
	def __init__(self,nodes,elementProperty):
		super(CBAR1,self).__init__('CBAR1',nodes,np.array([1.,0.,0.]),elementProperty);

		# Compute element rotation matrix T
		self.T = self.RotationMatrix()

	def __str__(self):
		return super(CBAR1,self).__str__()

	cpdef InitializeData(self,DataContainer dc):
		# Initialize stretch0 and stretch
		stretch0 = [1.0];
		stretch  = [1.0];

		self.stretch0ID =  dc.AddToData(dc.stretch0,dc.countStretch0,stretch0)
		self.stretchID  =  dc.AddToData(dc.stretch ,dc.countStretch ,stretch )
		return

	cpdef RotationMatrix(self):
		cdef np.ndarray[double, ndim=1] ex
		cdef np.ndarray[double, ndim=2] nodeT

		ex = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
		ex = ex/np.linalg.norm(ex)

		return spl.block_diag(ex,ex)

	cpdef BuildElementMatrix(self,DataContainer dc):
		cdef double l0,l, stretch
		cdef np.ndarray[double,ndim=1] ex0, ex
		cdef np.ndarray[double,ndim=2] mat
		cdef dict M

		# Compute stretch as ratio of direction vectors
		ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
		ex  = ex0 + np.array(self.nodes[1].Dof(dc))-np.array(self.nodes[0].Dof(dc))
		l0  = np.linalg.norm(ex0);
		l   = np.linalg.norm(ex);
		stretch = l/l0;
		stretch_rate = (stretch-eh.getStretch0(self,0,dc))/dc.dt;

		# Compute normal force and tangent
		M = self.property.Piola1Stiffness(stretch,stretch_rate,dc);

		# Material matrix in local coordinates
		mat = M['A']/l0*np.array([[1.,-1.],[-1.,1.]])

		# Transform back into global x,y coordinates
		mat = self.T.T.dot(mat).dot(self.T);

		return mat
		

	cpdef BuildInternalForceVector(self,DataContainer dc):
		cdef double l0,l, stretch
		cdef np.ndarray[double,ndim=1] ex0, ex, vec
		cdef dict M

		# Compute stretch as ratio of direction vectors
		ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
		ex  = ex0 + np.array(self.nodes[1].Dof(dc))-np.array(self.nodes[0].Dof(dc))
		l0  = np.linalg.norm(ex0);
		l   = np.linalg.norm(ex);
		stretch = l/l0;
		eh.setStretch(self,0,stretch,dc)
		stretch_rate = (stretch-eh.getStretch0(self,0,dc))/dc.dt;

		# Compute normal force and tangent
		M = self.property.Piola1Stiffness(stretch,stretch_rate,dc);

		# Residual force
		vec = M['P']*np.array([-1.,1.])

		 # Transform back into global x,y,z coordinates
		vec = self.T.T.dot(vec)

		return vec