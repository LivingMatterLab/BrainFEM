# -*- coding: utf-8 -*-
# Bar element for updated lagrangian method
from Element cimport *
from math import *
cimport ElementHelper as eh

cimport numpy as np
import numpy as np
import scipy.linalg as spl

cdef class CBEAM(Element):
	def __init__(self,nodes,elementProperty, restLength=0.0):
		super(CBEAM,self).__init__('CBEAM',nodes,np.array([1.,0.,0.]),elementProperty);

		self.nip = 1
		if(restLength==0.0):
			ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
			self.restLength  = np.linalg.norm(ex0);
		else:
			self.restLength = restLength

		# Compute element rotation matrix T
		self.T = self.RotationMatrix()

	def __str__(self):
		return super(CBEAM,self).__str__()

	cpdef InitializeData(self,DataContainer dc):
		# Initialize curvature and curvature0, theta and theta0
		curvLoc = [];
		for i in range(0,self.nip):
			curvLoc.append(np.zeros(3))

		self.curvature0ID =  dc.AddToData(dc.curvature0,dc.countCurvature0,curvLoc)
		self.curvatureID  =  dc.AddToData(dc.curvature ,dc.countCurvature ,curvLoc)
		self.theta0ID =  dc.AddToData(dc.theta0,dc.countTheta0,curvLoc)
		self.thetaID  =  dc.AddToData(dc.theta ,dc.countTheta ,curvLoc)
		

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

		ez = np.array([0.,0.,1.]);
		ez = ez - np.dot(ez,ex)*ex;
		if(np.linalg.norm(ez)<0.01):
			ez = np.array([1.,0.,0.]);
			ez = ez - np.dot(ez,ex)*ex;
		ez = ez/np.linalg.norm(ez);
		ey = np.cross(ez,ex);

		nodeT = np.squeeze([ex,ey,ez])

		return spl.block_diag(nodeT,nodeT,nodeT,nodeT)

	cpdef BuildElementMatrix(self,DataContainer dc):
		cdef int ndim, ip, i,j,p,q,r,s
		cdef double r_ip,s_ip,t_ip, detJ
		cdef np.ndarray[double,ndim=1] ve, veE,r0, r0E,dNx, wp,gp
		cdef np.ndarray[double,ndim=2]  mat,lam
		cdef dict deform, M

		udof  = np.array([0, 1, 2, 6,  7,  8]);	# disp dof
		wdof = np.array([3, 4, 5, 9, 10, 11]);	# curv dof

		ve  = self.Dof(dc);
		veE = self.T.dot(ve);        # ve in element coordinates
		r0  = self.Loc();
		r0E = self.T[np.ix_(udof,udof)].dot(r0);        # r0 in element coordinates

		if self.nip==1:
			gp = np.array([0.])
			wp = np.array([2.])
		elif self.nip==2:
			gp = 0.577350269189626* np.array([-1.,1.])
			wp = np.array([1.,1.]);

		
		mat = np.zeros((len(ve),len(ve)));
		ndim = 3;

		# Loop through integration points
		for ip in range(self.nip):
			# r coordinate of this integration point
			r_ip = gp[ip];

			# Compute deformation gradient
			deform = eh.LocalDeformation_B3(r_ip,r0E,veE[udof])
			N = deform['N']
			dNx = deform['dNx']
			detJ = deform['J'];

			# Compute tangent t1
			lam = spl.expm(eh.SkewToMat(eh.getTheta(self,ip,dc,ndim)))
			t1 = lam[0,:]

			# Compute tangent matrix and updated deformation gradient
			#M = self.property.Piola1Stiffness(deform['F']-t1,eh.getCurvature(self,ip,dc,ndim),dc.dt,ndim);
			M = self.property.CauchyStiffness(deform['F']-t1,eh.getCurvature(self,ip,dc,ndim),lam,dc.dt,ndim);

			# Compute necessary tangent matrices
			Cn = M['dPdF']
			Cm = M['dMdK']
			Fhat = eh.SkewToMat(deform['F'])
			CnF = Cn.dot(Fhat)
			FCnF = Fhat.dot(CnF)
			Phat = eh.SkewToMat(M['P'])
			Mhat = eh.SkewToMat(M['M'])
			FPhat = Fhat.dot(Phat)

			"""
			print 'Cn =  ' , Cn
			print 'CnF = ' , CnF
			print 'FCnF = ', FCnF
			print 'Cm  = ' , Cm

			print 'N   = ' , N 
			print 'dNx = ' , dNx
			"""
			

			# Loop through nodes and increment stiffness matrix
			for i in range(2):
				for j in range(2):
					eni = i*2*ndim;        # *2 due to additional rotations
					enj = j*2*ndim;        # *2 due to additional rotations
					for p in range(ndim):
						for q in range(ndim):
							# Material part
							mat[eni+p,     enj+q]      += dNx[i]*dNx[j]*Cn[p,q]*wp[ip]
							mat[eni+ndim+p,enj+q]      += N[i]*dNx[j]*CnF[p,q]*wp[ip]
							mat[eni+p,     enj+ndim+q] += dNx[i]*N[j]*CnF[q,p]*wp[ip]
							mat[eni+ndim+p,enj+ndim+q] += (dNx[i]*dNx[j]*Cm[p,q]-N[i]*N[j]*FCnF[p,q])*wp[ip]

							# Geometric part
							mat[eni+ndim+p,enj+q]      += N[i]*dNx[j]*Phat[q,p]*wp[ip]
							mat[eni+p,     enj+ndim+q] += dNx[i]*N[j]*Phat[p,q]*wp[ip]
							mat[eni+ndim+p,enj+ndim+q] += (N[i]*N[j]*FPhat[p,q]-dNx[i]*N[j]*Mhat[p,q])*wp[ip]
		
		# Multiply by detJ to correct for length
		mat*=detJ

		# Transform back into global x,y,z coordinates
		mat = self.T.T.dot(mat).dot(self.T);

		"""
		print 'mat[0,0] = ',mat[0,0]
		print 'mat[1,1] = ',mat[1,1]
		print 'mat[2,2] = ',mat[2,2]
		print 'mat[3,3] = ',mat[3,3]
		print 'mat[4,4] = ',mat[4,4]
		print 'mat[5,5] = ',mat[5,5]
		print 'mat[6,6] = ',mat[6,6]
		print 'mat[7,7] = ',mat[7,7]
		print 'mat[8,8] = ',mat[8,8]
		print 'mat[9,9] = ',mat[9,9]
		print 'mat[10,10] = ',mat[10,10]
		print 'mat[11,11] = ',mat[11,11]
		
		print 'mat[1,5]  = ',mat[1,5]
		print 'mat[7,11] = ',mat[7,11]
		print 'mat[2,4] = ',mat[2,4]
		print 'mat[8,10] = ',mat[8,10]
		print 'mat[4,5] = ',mat[4,5]
		print 'mat[10,11] = ',mat[10,11]
		print 'mat[1,7]  = ',mat[1,7]
		print 'mat[5,7]  = ',mat[5,7]
		print 'mat[1,11]  = ',mat[1,11]
		print 'mat[5,11]  = ',mat[5,11]
		"""

		return np.asarray(mat)
		
	cpdef BuildInternalForceVector(self,DataContainer dc):
		cdef int ndim, ip, i,j,p,q,r,s
		cdef double r_ip,s_ip,t_ip, detJ
		cdef np.ndarray[double,ndim=1] ve, veE, ve0, veE0,r0, r0E,dNx, wp,gp,vec,t1
		cdef np.ndarray[double,ndim=2] lam
		cdef dict deform, M

		udof  = np.array([0, 1, 2, 6,  7,  8]);	# disp dof
		wdof = np.array([3, 4, 5, 9, 10, 11]);	# curv dof

		ve  = self.Dof(dc);
		veE = self.T.dot(ve);        # ve in element coordinates
		ve0  = self.Dof0(dc);
		ve0E = self.T.dot(ve0);        # ve0 in element coordinates
		r0  = self.Loc();
		r0E = self.T[np.ix_(udof,udof)].dot(r0);        # r0 in element coordinates

		if self.nip==1:
			gp = np.array([0.])
			wp = np.array([2.])
		elif self.nip==2:
			gp = 0.577350269189626* np.array([-1.,1.])
			wp = np.array([1.,1.]);
		vec = np.zeros(len(ve));
		ndim = 3;

		# Compute increment in w at both nodes
		DW = np.zeros((ndim,2));
		# Node 0
		DW[:,0] = veE[3:6]-ve0E[3:6]
		# Node 1
		DW[:,1] = veE[9:12]-ve0E[9:12]

		# Loop through integration points
		for ip in range(self.nip):
			# r coordinate of this integration point
			r_ip = gp[ip];

			# Compute deformation gradient
			deform = eh.LocalDeformation_B3(r_ip,r0E,veE[udof])
			N = deform['N']
			dNx = deform['dNx']
			detJ = deform['J'];

            # Update curvature on this integration point
			DK = np.zeros(ndim);
			DWip = np.zeros(ndim);
			for j in range(2):
				DK += dNx[j]*DW[:,j];
				DWip += N[j]*DW[:,j];
			eh.setCurvature(self,ip,eh.getCurvature0(self,ip,dc,ndim)+DK,dc,ndim)
			self.UpdateTheta(ip,DWip,dc,ndim)

			# Compute tangent t1
			lam = spl.expm(eh.SkewToMat(eh.getTheta(self,ip,dc,ndim)))
			t1 = lam[0,:]
			
			"""
			#print 'lam = \n',lam
			print 'F  = ',deform['F']
			print 'lam = \n',lam
			print 't1 = ',t1
			print 'norm(t1) = ' , np.linalg.norm(t1)
			print 'F-t1 = ', deform['F']-t1
			print 'lam.t*F = ' , lam.dot(deform['F']-t1)
			"""

			# Compute tangent matrix and updated deformation gradient
			M = self.property.CauchyStiffness(deform['F']-t1,eh.getCurvature(self,ip,dc,ndim),lam,dc.dt,ndim);

			# Compute necessary vectors
			FP = np.cross(deform['F'],M['P']);

			# Loop through nodes and increment residual vector
			for i in range(2):
				eni = range(i*2*ndim,(i+1)*2*ndim); # *2 due to additional rotations
				vec[eni] += np.concatenate((dNx[i]*M['P'],\
					                        dNx[i]*M['M']-N[i]*FP))*wp[ip]
			"""
			print 'DK = ',DK
			print 'F  = ',deform['F']
			print 'P  = ',M['P']
			print 'FP = ',FP
			print 'Cn =  \n' , M['dPdF']
			print 'Cm =  \n' , M['dMdK']
			"""	
	
		# Multiply by detJ to correct for length
		vec*=detJ
		
		# Transform back into global x,y,z coordinates
		vec = self.T.T.dot(vec)

		"""
		print 'vec[0] = ', vec[0]
		print 'vec[1] = ', vec[1]
		print 'vec[2] = ', vec[2]
		print 'vec[3] = ', vec[3]
		print 'vec[4] = ', vec[4]
		print 'vec[5] = ', vec[5]
		print 'vec[6] = ', vec[6]
		"""

		return vec

	cpdef T_thw(self,np.ndarray[double,ndim=1] th_vec,int ndim):
		cdef np.ndarray[double,ndim=2] th_mat 
		cdef double th_norm, th_frac

		th_norm = np.linalg.norm(th_vec);
		
		if th_norm>1.e-6:
			th_frac = th_norm/2./np.tan(th_norm/2.);
			th_mat = eh.SkewToMat(th_vec)
			return th_frac*np.eye(ndim) - th_mat/2. + \
		           (1-th_frac)/th_norm/th_norm*np.outer(th_vec,th_vec)
		else:
			return np.eye(ndim)

	cpdef T_wth(self,np.ndarray[double,ndim=1] th_vec,int ndim):
		cdef np.ndarray[double,ndim=2] th_mat 
		cdef double th_norm, th_frac
		th_norm = np.linalg.norm(th_vec);

		if th_norm>1.e-6:
			th_frac = np.sin(th_norm)/th_norm;
			th_mat = eh.SkewToMat(th_vec)
			return th_frac*np.eye(ndim) + \
			       (1-np.cos(th_norm))/th_norm/th_norm*th_mat + \
		           (1-th_frac)/th_norm/th_norm*np.outer(th_vec,th_vec)
		else:
			return np.eye(ndim)

	cpdef UpdateTheta(self,int ip,np.ndarray[double,ndim=1] dw, DataContainer dc,int ndim):
		cdef np.ndarray[double,ndim=1] th0, th1, th 
		cdef int count

		# Get theta0
		th0 = eh.getTheta0(self,ip,dc,ndim)

		# Update theta in Newton raphson iteration
		th1 = th0
		th  = th0 + self.T_thw(th0,ndim).dot(dw)
		count = 0;
		while np.linalg.norm(th-th1)>1.e-6 and count<100:
			count+=1
			th1 = th 
			th  = th0 + self.T_thw(th1,ndim).dot(dw)
		
		# Check if converged, and update in that case.
		if count==100:
			raise "Theta in CBEAM.UpdateTheta did not converge after 100 iterations"
		else:
			# Update theta
			eh.setTheta(self,ip,th,dc,ndim)
			
