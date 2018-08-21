# -*- coding: utf-8 -*-
from Material cimport *
import numpy as np
cimport numpy as np
from math import *

cdef class MAT_NEOH(Material):
	def __init__(self,matProp):
		super(MAT_NEOH,self).__init__('MAT_NEOH');
		
		self.E      =   matProp['E'];
		self.nu     =   matProp['nu'];
		self.mu     = self.E/2.0/(1+self.nu);
		self.lamb   = self.E*self.nu / (1.0+self.nu)/(1.0-2.0*self.nu)
		try:
			self.D = matProp['D'];
		except:
			self.D = 0.

	def __str__(self):
		return super(MAT_NEOH,self).__str__() + "\t E = " + str(self.E)+"\t nu = " + str(self.nu)+ "\t D = " + str(self.D)

	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] Fe, int ndim):
		cdef double Je
		cdef np.ndarray[double, ndim=2] delta, Feinv, Pe
		cdef np.ndarray[double, ndim=4] Ae

		delta = np.eye(ndim);
		Je = np.linalg.det(Fe);
		Feinv = np.linalg.inv(Fe);

		# Error control for Je
		if(Je<=0.):
			Je = 0.001

		# 1st Piola Kirchhoff stress
		Pe = self.mu*Fe + (self.lamb*log(Je)-self.mu)*np.transpose(Feinv);

		# 4th order elastic stiffness Ae = dPe/dFe
		Ae = np.empty((ndim,ndim,ndim,ndim))
		for i in range(ndim):
			for j in range(ndim):
				for k in range(ndim):
					for l in range(ndim):
						Ae[i,j,k,l] = self.lamb                     * Feinv[j,i]*Feinv[l,k] \
							- (self.lamb*log(Je)-self.mu) * Feinv[l,i]*Feinv[j,k] \
							+ self.mu                     * delta[i,k]*delta[j,l];

		return {'Pe':Pe,'Ae':Ae}

	cpdef Piola1Stiffness_1d(self,double stretche, double stretche_rate = 0.0, double dt = 0.0):
		# 1st Piola Kirchhoff stress
		Pe = self.mu*(stretche-1./(stretche**2))

		# Elastic stiffness Ae = dPe/dstretche
		Ae = self.mu*(1.+2./(stretche**3))
		
		return {'Pe':Pe,'Ae':Ae}

	cpdef CauchyStiffness_1d(self, double stretche, double stretche_rate = 0.0, double dt = 0.0):
		# Cauchy stress
		sige = self.mu*(stretche**2-1./stretche)

		# Elastic stiffness Ce = dsige/dstretche
		Ce = self.mu*(2.*stretche+1./(stretche**2))
		
		return {'sige':sige,'Ce':Ce}