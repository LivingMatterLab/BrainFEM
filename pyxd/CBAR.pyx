# -*- coding: utf-8 -*-
# Bar element for updated lagrangian method
from Element cimport *
from math import *
cimport ElementHelper as eh

cimport numpy as np
import numpy as np
import scipy.linalg as spl

cdef class CBAR(Element):
	def __init__(self,nodes,elementProperty, restLength=0.0):
		super(CBAR,self).__init__('CBAR',nodes,np.array([1.,0.,0.]),elementProperty);

		if(restLength==0.0):
			ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
			self.restLength  = np.linalg.norm(ex0);
		else:
			self.restLength = restLength

	def __str__(self):
		return super(CBAR,self).__str__()

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
		cdef double l, stretch, stretch_rate
		cdef np.ndarray[double,ndim=1] ex0, ex
		cdef np.ndarray[double,ndim=2] nn,Inn,mat,T,mat2
		cdef int nZeroPad
		cdef dict M

		# Compute stretch as ratio of direction vectors
		ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
		ex  = ex0 + np.array(self.nodes[1].DispDof(dc))-np.array(self.nodes[0].DispDof(dc))
		l   = np.linalg.norm(ex);
		stretch = l/self.restLength;
		stretch_rate = (stretch-eh.getStretch0(self,0,dc))/dc.dt;

		"""
		# Compute normal force and tangent
		M = self.property.CauchyStiffness(stretch,stretch_rate,dc.dt);

		# Material matrix in local coordinates
		mat = ( M['C']/l0 - M['sig']/l)*np.array([[1.,-1.],[-1.,1.]])
		"""
		# Compute normal force and tangent
		M = self.property.Piola1Stiffness(stretch,stretch_rate,dc);

		
		# Compute element matrix
		ex = ex/l
		nn = np.outer(ex,ex)
		Inn = np.eye(len(ex))-nn

		# Check if we need to pad zeros to matrix (e.g. if the 
		# node also contains rotation dof)
		nZeroPad = (self.datRID[1]-self.datRID[0]-2*len(nn))/2
		if nZeroPad>0:
			nn  = spl.block_diag(nn, np.zeros((nZeroPad,nZeroPad)))
			Inn = spl.block_diag(Inn,np.zeros((nZeroPad,nZeroPad)))

		mat = M['A']/self.restLength * np.bmat([[ nn, -1* nn],[-1* nn, nn]]) + \
		      M['P']/l  * np.bmat([[Inn, -1*Inn],[-1*Inn,Inn]])
		
		return np.asarray(mat)
		
	cpdef BuildInternalForceVector(self,DataContainer dc):
		cdef double l, stretch, stretch_rate
		cdef np.ndarray[double,ndim=1] ex0, ex, vec
		cdef int nZeroPad
		cdef dict M

		# Compute stretch as ratio of direction vectors
		ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
		ex  = ex0 + np.array(self.nodes[1].DispDof(dc))-np.array(self.nodes[0].DispDof(dc))
		l   = np.linalg.norm(ex);
		stretch = l/self.restLength;
		eh.setStretch(self,0,stretch,dc)
		stretch_rate = (stretch-eh.getStretch0(self,0,dc))/dc.dt;

		"""
		# Compute normal force and tangent
		M = self.property.CauchyStiffness(stretch,stretch_rate,dc.dt);

		# Residual force
		vec = M['sig']*np.array([-1.,1.])
		"""

		# Compute normal force and tangent
		M = self.property.Piola1Stiffness(stretch,stretch_rate,dc);

		# Check whether we need to pad zeros to matrix (e.g. if the 
		# node also contains rotation dof)
		nZeroPad = (self.datRID[1]-self.datRID[0]-2*len(ex))/2
		if nZeroPad>0:
			ex = np.concatenate((ex,np.zeros(nZeroPad)))

		# Compute residual force vector
		ex = ex/l
		vec = M['P']* np.concatenate((-1*ex,ex))

		return vec