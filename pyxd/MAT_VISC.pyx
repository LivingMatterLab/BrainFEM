# -*- coding: utf-8 -*-
from Material cimport *
import numpy as np
cimport numpy as np
from math import *

cdef class MAT_VISC(Material):
	def __init__(self,matProp):
		super(MAT_VISC,self).__init__('MAT_VISC');
		
		self.eta      =   matProp['eta'];

	def __str__(self):
		return super(MAT_VISC,self).__str__() + "\t eta = " + str(self.eta)

	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] Fe, int ndim):
		raise "Piola1Stiffness in MAT_VISC not implemented yet"

	cpdef Piola1Stiffness_1d(self,double stretche, double stretche_rate, double dt):
		# 1st Piola Kirchhoff stress
		Pe = self.eta*stretche_rate
		
		# Elastic stiffness Ae = dPe/dstretche
		Ae = self.eta/dt
		
		return {'Pe':Pe,'Ae':Ae}

	cpdef CauchyStiffness_1d(self, double stretche, double stretche_rate, double dt):
		# Cauchy stress
		sige = self.eta*stretche_rate

		# Elastic stiffness Ce = dsige/dstretche
		Ce = self.eta/dt
		
		return {'sige':sige,'Ce':Ce}
	"""
	cpdef Piola1Stiffness_1d(self,double stretche, double stretche_rate, double dt):
		# 1st Piola Kirchhoff stress
		Pe = 2.*self.eta*stretche_rate/(stretche**2)

		# Elastic stiffness Ae = dPe/dstretche
		Ae = -2.*Pe/stretche+2.*self.eta/(stretche**2)/dt
		
		return {'Pe':Pe,'Ae':Ae}

	cpdef CauchyStiffness_1d(self, double stretche, double stretche_rate, double dt):
		# Cauchy stress
		sige = 2.*self.eta*stretche_rate/stretche

		# Elastic stiffness Ce = dsige/dstretche
		Ce = -sige/stretche+2.*self.eta/stretche/dt
		
		return {'sige':sige,'Ce':Ce}
	"""