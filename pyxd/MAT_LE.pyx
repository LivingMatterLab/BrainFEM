# -*- coding: utf-8 -*-
from Material cimport *
import numpy as np
cimport numpy as np
from math import *

cdef class MAT_LE(Material):
	def __init__(self,matProp):
		super(MAT_LE,self).__init__('MAT_LE');
		
		self.E      =   matProp['E'];
		self.nu     =   matProp['nu'];
		self.mu     = self.E/2.0/(1+self.nu);
		self.lamb   = self.E*self.nu / (1.0+self.nu)/(1.0-2.0*self.nu)
		try:
			self.D = matProp['D'];
		except:
			self.D = 0.

		print self.mu

	def __str__(self):
		return super(MAT_LE,self).__str__() + "\t E = " + str(self.E)+"\t nu = " + str(self.nu)

	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] Fe, int ndim):
		cdef double Je
		cdef np.ndarray[double, ndim=2] beta, eps, delta, Pe
		cdef np.ndarray[double, ndim=4] Ae

		delta = np.eye(ndim);

		# Compute small strain
		beta = Fe-delta;				# displacement gradient
		eps = 0.5*(beta.T+beta)			# small strain tensor

		# 1st Piola Kirchhoff stress (or any other stress measure)
		Pe = self.lamb*np.trace(eps)*delta + 2.*self.mu*eps
		
		# 4th order elastic stiffness Ae = dPe/dFe
		Ae = np.empty((ndim,ndim,ndim,ndim))
		for i in range(ndim):
			for j in range(ndim):
				for k in range(ndim):
					for l in range(ndim):
						Ae[i,j,k,l] = self.lamb * delta[i,j]*delta[k,l] \
									+ self.mu   * (delta[i,k]*delta[j,l]+delta[i,l]*delta[j,k]);
		return {'Pe':Pe,'Ae':Ae}

	cpdef Piola1Stiffness_1d(self,double stretche, double stretche_rate = 0.0, double dt = 0.0):
		# 1st Piola Kirchhoff stress
		Pe = self.E*(stretche-1.0)

		# Elastic stiffness Ae = dPe/dstretche
		Ae = self.E
		
		return {'Pe':Pe,'Ae':Ae}

	cpdef CauchyStiffness_1d(self, double stretche, double stretche_rate = 0.0, double dt = 0.0):
		# Cauchy stress
		sige = self.E*(stretche-1.0)

		# Elastic stiffness Ce = dsige/dstretche
		Ce = self.E
		
		return {'sige':sige,'Ce':Ce}