# -*- coding: utf-8 -*-
from Material cimport *
import numpy as np
import scipy.linalg as spl
cimport numpy as np
from math import *

cdef class MAT_OGDEN(Material):
	def __init__(self,matProp):
		super(MAT_OGDEN,self).__init__('MAT_OGDEN');
		
		self.mu_i      = matProp['mu']
		self.alpha_i   = matProp['alpha']
		self.nu     =   matProp['nu'];

		self.mu = 0;
		for i in range(len(self.mu_i)):
			self.mu+=0.5*self.mu_i[i]*self.alpha_i[i]

		self.E      = self.mu*2.0*(1+self.nu);
		#self.K      = self.E/3.0/(1-2*self.nu);
		self.lamb   = self.E*self.nu / (1.0+self.nu)/(1.0-2.0*self.nu)
		self.beta   = self.nu/(1.-2.*self.nu) 
		
		try:
			self.D = matProp['D'];
		except:
			self.D = 0.

	def __str__(self):
		return super(MAT_OGDEN,self).__str__() +"\t mu_i = [" + ', '.join(str(mu) for mu in self.mu_i) + "]"+ "\t alpha_i = [" + ', '.join(str(a) for a in self.alpha_i) + "]\t nu = " + str(self.nu)+ "\t D = " + str(self.D)

	cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] Fe, int ndim):

		cdef int eid, oid, nd
		cdef double Je,

		cdef np.ndarray[double, ndim=1] eigVal, omega, domega, ddomega, aa, bb, IC
		cdef np.ndarray[double, ndim=2] eigVec, cc, Feinv, Pe,delta, Ce
		cdef np.ndarray[double, ndim=4] Ae
		cdef list dIC, ddIC

		# Compute inverse and eigenvalues
		Je = np.linalg.det(Fe);
		nd = 3;	
		delta = np.eye(nd)
		if ndim==2:
			Fe     = spl.block_diag(Fe,np.array([1.]))
		Feinv = np.linalg.inv(Fe)
		Ce    = Fe.T.dot(Fe)
		eigVal, eigVec = np.linalg.eigh(Ce)

		# Invariants
		IC = np.zeros(nd)
		IC[0] = np.trace(Ce)
		IC[1] = 0.5*(np.trace(Ce)**2-np.trace(Ce.dot(Ce)))
		IC[2] = np.linalg.det(Ce)

		# Compute omega and its derivatives
		omega   = np.zeros(nd);			#omega_1, omega_2, omega_3
		domega  = np.zeros(nd+1);		#domega_1/dlambda_1 .. d(omega_1+omega_2+omega_3)/dJ
		ddomega = np.zeros(nd+1);       #d2omega_1/d2lambda_1 .. d2(omega_1+omega_2+omega_3)/d2J

		for eid in range(nd):
			for oid in range(len(self.mu_i)):


				omega[eid]+= self.mu_i[oid]/self.alpha_i[oid]* \
					         (eigVal[eid]**(self.alpha_i[oid]/2.)-1.)


				domega[eid]+= self.mu_i[oid]/2.* \
					          (eigVal[eid])**(self.alpha_i[oid]/2.-1.)


				ddomega[eid]+= self.mu_i[oid]/2.*(self.alpha_i[oid]/2.-1)* \
					          eigVal[eid]**(self.alpha_i[oid]/2.-2.)

			    
				domega[nd]  -= self.mu_i[oid]/nd* \
					         Je**(-self.beta*self.alpha_i[oid]-1.)

				ddomega[nd] += self.mu_i[oid]/nd*(self.beta*self.alpha_i[oid]+1.) * \
					         Je**(-self.beta*self.alpha_i[oid]-2.)

		
		aa,bb,cc = self.GetCoefficients(eigVal,IC,omega,domega,ddomega,nd)

		# Compute derivatives of invariants
		dIC = [];
		dIC.append(2*Fe);						# = dI1/dFe
		dIC.append(2*IC[0]*Fe-2*Fe.dot(Ce))		# = dI2/dFe
		dIC.append(2*IC[2]*Feinv.T)  			# = dI3/dFe
		dIC.append(Je*Feinv.T)  				# = dJ/dFe

		ddIC = [];
		ddIC.append(2*np.einsum('ik,jl->ijkl',delta,delta))

		ddIC.append(2*IC[0]*np.einsum('ik,jl->ijkl',delta,delta) \
		             +4*np.einsum('ij,kl->ijkl',Fe,Fe) \
		             -2*np.einsum('ik,jl->ijkl',delta,Ce) \
		             -2*np.einsum('il,jk->ijkl',Fe,Fe) \
		             -2*np.einsum('ik,jl->ijkl',Fe.dot(Fe.T),delta))

		ddIC.append(2*IC[2]*(2*np.einsum('lk,ji->ijkl',Feinv,Feinv) \
			                  -np.einsum('li,jk->ijkl',Feinv,Feinv)))


		ddIC.append(     Je*(  np.einsum('lk,ji->ijkl',Feinv,Feinv) \
			                  -np.einsum('li,jk->ijkl',Feinv,Feinv)))


		# 1st Piola Kirchhoff stress
		Pe = np.zeros((nd,nd))#0.5*self.K*(-Je**(-2)+1.)*Feinv.T		# Compressible
		for eid in range(nd):
			Pe +=aa[eid]*dIC[eid]
		Pe+=domega[nd]*dIC[nd]		#domega/dJ*dJ/dF

		# 4th order elastic stiffness Ae = dPe/dFe
		#Ae = self.K/(Je**2)*np.einsum('lk,ji->ijkl',Feinv,Feinv) \
		#	 -0.5*self.K*(-Je**(-2.)+1)*np.einsum('li,jk->ijkl',Feinv,Feinv)	# Compressible
		Ae = np.zeros((nd,nd,nd,nd))	

		for eid in range(nd):
			Ae +=bb[eid]*ddIC[eid]
			for oid in range(nd):
				Ae +=cc[eid,oid]*np.einsum('ij,kl->ijkl',dIC[eid],dIC[oid])

		Ae+=domega[nd]*ddIC[nd]		# domega/dJ*d2J/dF2
		Ae+=ddomega[nd]*np.einsum('ij,kl->ijkl',dIC[nd],dIC[nd])	# d2/omage/dJ2*(dJ/dF)**2


		if ndim==2:
			Pe = Pe[np.ix_([0,1],[0,1])]
			Ae = Ae[np.ix_([0,1],[0,1],[0,1],[0,1])]

		return {'Pe':Pe,'Ae':Ae}

	cpdef Piola1Stiffness_1d(self,double stretche, double stretche_rate = 0.0, double dt = 0.0):
		raise 'ERROR: Piola1Stiffness_1d not implemented for Gent material'

	cpdef CauchyStiffness_1d(self, double stretche, double stretche_rate = 0.0, double dt = 0.0):
		raise 'ERROR: Piola1Stiffness_1d not implemented for Gent material'


	cpdef GetCoefficients(self, np.ndarray[double,ndim=1] eigVal, \
	                            np.ndarray[double,ndim=1] IC,\
	                            np.ndarray[double,ndim=1] omega, \
	                            np.ndarray[double,ndim=1] domega, \
	                            np.ndarray[double,ndim=1] ddomega, \
	                            int ndim):
		
		cdef np.ndarray[double, ndim=1] aa, bb,dd
		cdef np.ndarray[double, ndim=2]  cc
		cdef list duplicity
		

		# Compute duplicity
		duplicity = [];
		if np.abs(eigVal[1]-eigVal[0])<1.e-6:
			duplicity.append(0);
		if np.abs(eigVal[2]-eigVal[1])<1.e-6:
			duplicity.append(1);
		#print eigVal
		#print duplicity

		# Compute aa, bb, cc coefficients
		aa   = np.zeros(ndim);
		bb   = np.zeros(ndim);
		cc   = np.zeros((ndim,ndim));
		
		if len(duplicity)==0:		# No duplicate eigenvalues
			dd = np.zeros(ndim)
			dd[0] = (eigVal[0]-eigVal[1])*(eigVal[0]-eigVal[2])
			dd[1] = (eigVal[1]-eigVal[2])*(eigVal[1]-eigVal[0])
			dd[2] = (eigVal[2]-eigVal[0])*(eigVal[2]-eigVal[1])
			for eid in range(ndim):
				aa[0]  += domega[eid]*eigVal[eid]**2/dd[eid]
				aa[1]  -= domega[eid]*eigVal[eid]   /dd[eid]
				aa[2]  += domega[eid]               /dd[eid]

				cc[0,0]+= eigVal[eid]**2/dd[eid]**2 * \
				          (ddomega[eid]*eigVal[eid]**2 + \
				           2*domega[eid]*(3*IC[2]-eigVal[eid]*IC[1])/dd[eid])


				cc[1,1]+= 1./dd[eid]**2 * \
				          (ddomega[eid]*eigVal[eid]**2 + \
				           2*domega[eid]*(IC[2]-eigVal[eid]**3)/dd[eid])


				cc[2,2]+= 1./dd[eid]**2 * \
				          (ddomega[eid] + \
				           2*domega[eid]*(IC[0]-3.*eigVal[eid])/dd[eid])

				cc[0,1]+= eigVal[eid]/dd[eid]**2 * \
				          (-ddomega[eid]*eigVal[eid]**2 + \
				           domega[eid]*(eigVal[eid]**2*IC[0]-3*IC[2])/dd[eid])

				cc[1,2]+= 1./dd[eid]**2 * \
				          (-ddomega[eid]*eigVal[eid] + \
				           domega[eid]*(3*eigVal[eid]**2-IC[1])/dd[eid])

			bb[0]=aa[0]; bb[1]=aa[1]; bb[2]=aa[2];
			cc[0,2] = cc[1,1]
		elif len(duplicity)==1:		# One duplicate eigenvalues
			# We want eigVal[1]==eigVal[2] So switch order if not the case
			if duplicity[0]==0:
				print 'Warning: Check implementation of reordering in MAT_OGDEN!!!!'
				eigVal = eigVal[np.ix_([2,1,0])]

			dd = np.zeros(2)
			dd[0] = domega[0]-domega[1]-(eigVal[0]-eigVal[1])*ddomega[1]
			dd[1] = ddomega[0]-ddomega[1]

			aa[2] = (domega[0]-domega[1])/(eigVal[0]-eigVal[1])**2
			aa[1] =-aa[2]*eigVal[0]
			aa[0] =-aa[1]*eigVal[0]+domega[1]

			bb[2] = dd[0]/(eigVal[0]-eigVal[1])**2
			bb[1] =-bb[2]*eigVal[0]-ddomega[1]
			bb[0] = bb[2]*eigVal[0]**2+domega[1]+ddomega[1]*(eigVal[0]+eigVal[1])

			cc[0,0] = -omega[0]*4*eigVal[0]**3*eigVal[1]/(eigVal[0]-eigVal[1])**5 +\
			           omega[1]*(eigVal[0]/(eigVal[0]-eigVal[1]))**4 + \
			           ddomega[1]

			cc[1,1] = -omega[0]*2*eigVal[0]*(eigVal[0]+eigVal[1])/(eigVal[0]-eigVal[1])**5 +\
			           omega[1]*eigVal[0]**2/(eigVal[0]-eigVal[1])**4

			cc[2,2] = -omega[0]*4/(eigVal[0]-eigVal[1])**5 +\
			           omega[1]/(eigVal[0]-eigVal[1])**4

			cc[0,1] =  omega[0]*eigVal[0]**2*(eigVal[0]+3*eigVal[1])/(eigVal[0]-eigVal[1])**5 +\
			          -omega[1]*eigVal[0]**3/(eigVal[0]-eigVal[1])**4

			cc[0,2] =  cc[1,1]

			cc[1,2] =  omega[0]*(3*eigVal[0]+eigVal[1])/(eigVal[0]-eigVal[1])**5 +\
			          -omega[1]*eigVal[0]/(eigVal[0]-eigVal[1])**4

		else:		# All eigenvalues the same
			aa[0] = domega[0]
			bb[0] = domega[0]+2*eigVal[0]*ddomega[0]
			bb[1] = -ddomega[0]
			cc[0,0] = ddomega[0]

		if True:
			cc[1,0] = cc[0,1]
			cc[2,1] = cc[1,2]
			cc[2,0] = cc[0,2]

		return (aa,bb,cc)