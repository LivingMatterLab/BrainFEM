# -*- coding: utf-8 -*-
import numpy as np
cimport numpy as np
import scipy.optimize as sco

cdef class CONTACT:
	def __init__(self, slaveNodes, masterNodes,penaltyParameter):
		self.slaveNodes     = slaveNodes
		self.masterNodes    = masterNodes
		self.numSlaveNodes  = len(self.slaveNodes)
		self.numMasterNodes = len(self.masterNodes)
		self.penaltyParameter = penaltyParameter


	def __str__(self):
		return "CONTACT: \t SlaveNodes = ["+ ', '.join(str(x.localID) for x in self.slaveNodes)+"]" \
		               +"\t MasterNodes = ["+ ', '.join(str(x.localID) for x in self.masterNodes)+"]"

	cpdef BuildMatrix(self, DataContainer dc):
		ndim = 2
		mat = np.zeros((self.datRID[1]-self.datRID[0],self.datRID[1]-self.datRID[0]));

		# Collect current and old locations of all master nodes
		mLoc0 = [];
		mLoc = [];
		for m in self.masterNodes:
			mLoc0.append(np.array(m.loc)+np.array(m.DispDof0(dc)))
			mLoc.append(np.array(m.loc)+np.array(m.DispDof(dc)))

		# Compute current normal vectors to each master line element
		tangVec = [];normVec = [];lenElem = []
		for im in range(self.numMasterNodes-1):
			ex = mLoc[im+1]-mLoc[im]
			lenElem.append(np.linalg.norm(ex))
			ex = ex/lenElem[-1]
			tangVec.append(ex)
			normVec.append(np.array([-ex[1],ex[0]]))

		# Loop through slave nodes and check for contact with master elements
		vid_S = -ndim;
		for sID in range(self.numSlaveNodes):
			s = self.slaveNodes[sID]
			sLoc  = np.array(s.loc)+np.array(s.DispDof(dc))
			sLoc0 = np.array(s.loc)+np.array(s.DispDof0(dc))
			vid_S+=ndim

			# Check contact of slave node
			contactInfo = self.checkSlaveContact(sLoc,mLoc,lenElem,tangVec,normVec)
			penDepth = contactInfo['penDepth']
			penID    = contactInfo['penID']
			penXi    = contactInfo['penXi']
			penNVec  = contactInfo['penNVec']
			penTVec  = contactInfo['penTVec']
			penElLen = contactInfo['penElLen']

			# Check if penetration occurred and act accordingly
			if penDepth>0:
				vid_M = (self.numSlaveNodes+penID)*ndim
				penForce, dpenForce = self.getPenalty(penDepth)


				# Compute normal/tangent vectors
				Ns  = np.concatenate([penNVec,-(1.-penXi)*penNVec,-penXi*penNVec])
				N0s = np.concatenate([np.zeros(ndim),           -penNVec,       penNVec])

				Ts  = np.concatenate([penTVec,-(1.-penXi)*penTVec,-penXi*penTVec])
				T0s = np.concatenate([np.zeros(ndim),           -penTVec,      penTVec])
				
				# Add penForce to vec
				kid = [i for i in range(vid_S,vid_S+ndim)]
				kid.extend([i for i in range(vid_M,vid_M+2*ndim)])

				Kd = -(np.outer(N0s,Ts)+np.outer(Ts,N0s) +\
					   penDepth/penElLen*np.outer(N0s,N0s))/penElLen

				"""
				print 'Kd=\n',penForce*Kd
				print 'NN=\n',2*self.penaltyParameter*penDepth*np.outer(Ns,Ns)
				print '\n'
				"""

				mat[np.ix_(kid,kid)]+= dpenForce*np.outer(Ns,Ns)#+penForce*Kd

				# Add some stability (not consistent)
				#mat[np.ix_(kid,kid)]+= 1.e-5*dpenForce*np.outer(Ts,Ts)
			
		return mat

	cpdef BuildVector(self, DataContainer dc):
		cdef np.ndarray[double,ndim=1] loadVec

		ndim = 2;
		vec = np.zeros(self.datRID[1]-self.datRID[0]);

		# Collect current and old locations of all master nodes
		mLoc0 = [];
		mLoc = [];
		for m in self.masterNodes:
			mLoc0.append(np.array(m.loc)+np.array(m.DispDof0(dc)))
			mLoc.append(np.array(m.loc)+np.array(m.DispDof(dc)))

		# Compute current normal vectors to each master line element
		tangVec = [];normVec = [];lenElem = []
		for im in range(self.numMasterNodes-1):
			ex = mLoc[im+1]-mLoc[im]
			lenElem.append(np.linalg.norm(ex))
			ex = ex/lenElem[-1]
			tangVec.append(ex)
			normVec.append(np.array([-ex[1],ex[0]]))

		# Loop through slave nodes and check for contact with master elements
		vid_S = -ndim;
		for sID in range(self.numSlaveNodes):
			s = self.slaveNodes[sID]
			sLoc  = np.array(s.loc)+np.array(s.DispDof(dc))
			sLoc0 = np.array(s.loc)+np.array(s.DispDof0(dc))
			vid_S+=ndim

			# Check contact of slave node
			contactInfo = self.checkSlaveContact(sLoc,mLoc,lenElem,tangVec,normVec)
			penDepth = contactInfo['penDepth']
			penID    = contactInfo['penID']
			penXi    = contactInfo['penXi']
			penNVec  = contactInfo['penNVec']
			penTVec  = contactInfo['penTVec']
			penElLen = contactInfo['penElLen']

			"""
			if s.localID==38:
				print contactFound
				print penDepth
				print penID
				print penXi
				print penVec
				print penElLen
				print '\n'
			"""

			# Check if penetration occurred and act accordingly
			if penDepth>0:
				#print s.localID
				vid_M = (self.numSlaveNodes+penID)*ndim
				penForce, dpenForce = self.getPenalty(penDepth)

				# Compute normal/tangent vectors
				Ns  = np.concatenate([penNVec,-(1.-penXi)*penNVec,-penXi*penNVec])
				
				# Add penForce to vec
				vid = [i for i in range(vid_S,vid_S+ndim)]
				vid.extend([i for i in range(vid_M,vid_M+2*ndim)])
				vec[vid]+=penForce*Ns

				"""
				print sLoc
				print mLoc[penID], mLoc[penID+1]
				print self.masterNodes[penID].DispDof(dc)
				print mLoc[penID]+self.masterNodes[penID].DispDof0(dc)
				print 'Slave node penetrated: ', penDepth
				print 'Slave node id:         ', vid_S
	
				print 'Master node id:        ', vid_M
				#print vec
				"""

		return vec

	cpdef getHermiteG(self, double g,double alpha,\
		           np.ndarray[double,ndim=1] x1,\
		           np.ndarray[double,ndim=1] x2,\
		           np.ndarray[double,ndim=1] x3,\
		           np.ndarray[double,ndim=1] x4):
		cdef np.ndarray[double,ndim=1] tg, xg 


		xg = 0.25*alpha* ( (g**2-2*g+1)*x2 + (g**2+2*g+1)*x4  ) + \
		     (1.-0.5*alpha*(g*g+1))*x3

		tg = alpha* ( (g-1.)/2.*x2 - g*x3 + (g+1.)/2.*x4)

		return np.inner((x1-xg),tg)

	cpdef getHermiteDepth(self, double g,double alpha,\
		           np.ndarray[double,ndim=1] x1,\
		           np.ndarray[double,ndim=1] x2,\
		           np.ndarray[double,ndim=1] x3,\
		           np.ndarray[double,ndim=1] x4):
		cdef np.ndarray[double,ndim=1] tg,ng, xg 


		xg = 0.25*alpha* ( (g**2-2*g+1)*x2 + (g**2+2*g+1)*x4  ) + \
		     (1.-0.5*alpha*(g*g+1))*x3

		tg = alpha* ( (g-1.)/2.*x2 - g*x3 + (g+1.)/2.*x4)
		tg = tg/np.linalg.norm(tg)
		ng = np.array([-tg[1],tg[0]])

		return (np.inner((x1-xg),ng), ng,tg)


	cpdef checkSlaveContact(self, np.ndarray[double,ndim=1] sLoc, list mLoc, list lenElem, list tangVec, list normVec):
		penDepth = -np.inf
		penID    = -1
		penXi    = 0.5
		penNVec  = np.array([0.,1.])
		penTVec  = np.array([1.,0.])
		penElLen = 1.
		im = 0;

		alpha = 0.1
		epsXi = alpha/2.

		for im in range(self.numMasterNodes-1):


			# Compute xi, indicating whether slave node is 'inside' two master nodes
			xi = np.inner(sLoc-mLoc[im],tangVec[im])/lenElem[im]

			if xi>-epsXi and xi<1+epsXi:
				# Compute vector from center element to slave node
				#ems  = sLoc-0.5*(mLoc[im] +mLoc[im+1])
				ems  = sLoc-(1-xi)*mLoc[im]-xi*mLoc[im+1]
				lms  = np.inner(normVec[im],ems)

				if lms < -penDepth:
					penDepth = -lms
					penID = im
					penXi = xi
					penNVec = normVec[penID]
					penTVec = tangVec[penID]
					penElLen = lenElem[penID]

		if penDepth < -1.e9:
			print 'WARNING: slave node not visible by any master element!!'
			penDepth = 0;

		# Check if Hermite interpolation is required
		if penXi>-epsXi and penXi<=alpha and not penID==0:
			# Hermite interpolation with previous element
			gg = sco.fsolve(self.getHermiteG, xi/alpha, args=(alpha,sLoc,mLoc[penID-1],mLoc[penID],mLoc[penID+1]) )[0]
			dep,nvec,tvec = self.getHermiteDepth(gg,alpha,sLoc,mLoc[penID-1],mLoc[penID],mLoc[penID+1])
			
			if gg<-1 or gg>1:
				print 'WARNING: gg = ', gg
			
			penDepth = -dep
			penElLen = 0.5*(gg+1)*lenElem[penID] + \
			           0.5*(1-gg)*lenElem[penID-1]
			
			penID    = penID-1         if gg<0 else penID
			penXi    = alpha*gg+1.  if gg<0 else alpha*gg
			penNVec   = nvec
			penTVec   = tvec
			

			#print 'normVec - = ', normVec[penID]
			#print 'penVec    = ', penNVec
			#print 'normVec + = ', normVec[penID+1]
			#print 'penDepth  = ', penDepth
			#print 'penXi     = ', penXi
			#print '\n'
				
		elif penXi>=1.-alpha and penXi<1.+epsXi and not penID==(self.numMasterNodes-2):
			# Hermite interpolation with next element
			
			gg = sco.fsolve(self.getHermiteG, (xi-1.)/alpha, args=(alpha,sLoc,mLoc[penID],mLoc[penID+1],mLoc[penID+2]) )[0]
			dep,nvec,tvec = self.getHermiteDepth(gg,alpha,sLoc,mLoc[penID],mLoc[penID+1],mLoc[penID+2])

			if gg<-1 or gg>1:
				print 'WARNING: gg = ', gg
			
			penDepth = -dep
			penElLen = 0.5*(gg+1)*lenElem[penID+1] + \
			           0.5*(1-gg)*lenElem[penID]
			penID    = penID           if gg<0 else penID+1
			penXi    = alpha*gg+1.  if gg<0 else alpha*gg
			penNVec   = nvec
			penTVec   = tvec

			#print 'normVec - = ', normVec[penID]
			#print 'penVec    = ', penNVec
			#print 'normVec + = ', normVec[penID+1]
			#print 'penDepth  = ', penDepth
			#print 'penXi     = ', penXi
			#print '\n'
			

	

		return {'penDepth':penDepth,'penID':penID,'penXi':penXi,'penNVec':penNVec,'penTVec':penTVec,'penElLen':penElLen}


	"""
	cpdef checkSlaveContact(self, np.ndarray[double,ndim=1] sLoc, list mLoc, list lenElem, list tangVec, list normVec):
		penDepth = np.inf
		penID    = -1
		penXi    = 0.5
		penNVec  = np.array([0.,0.])
		penTVec  = np.array([0.,0.])
		penElLen = 0.
		contactFound = False
		im = 0;

		while im<(self.numMasterNodes-1) and not contactFound:


			# Compute xi, indicating whether slave node is 'inside' two master nodes
			xi = np.inner(sLoc-mLoc[im],tangVec[im])/lenElem[im]

			# Compute vector from center element to slave node
			#ems  = sLoc-0.5*(mLoc[im] +mLoc[im+1])
			ems  = sLoc-(1-xi)*mLoc[im]-xi*mLoc[im+1]
			lms  = np.inner(normVec[im],ems)

			
			# Check if conditions for contact are met
			alpha = 0.1
			epsXi = alpha/2.
			if lms <0 and ( (xi>alpha and xi<1.-alpha) or \
				             (xi>0 and xi<alpha and im==0) or \
				             (xi>1.-alpha and xi<1 and im==self.numMasterNodes-2)):
				# Penetration in 'center' of element
				contactFound = True
				penDepth = -lms
				penID = im
				penXi = xi
				penNVec = normVec[penID]
				penTVec = tangVec[penID]
				penElLen = lenElem[penID]

			elif lms <0.2*lenElem[im] and xi>-epsXi and xi<=alpha and not im==0:
				# Hermite interpolation with previous element
				gg = sco.fsolve(self.getHermiteG, xi/alpha, args=(alpha,sLoc,mLoc[im-1],mLoc[im],mLoc[im+1]) )[0]
				dep,nvec,tvec = self.getHermiteDepth(gg,alpha,sLoc,mLoc[im-1],mLoc[im],mLoc[im+1])
				
				if gg<-1 or gg>1:
					print 'WARNING: gg = ', gg

				if dep < 0:
					contactFound = True
					penDepth = -dep
					penID    = im-1         if gg<0 else im
					penXi    = alpha*gg+1.  if gg<0 else alpha*gg
					penNVec   = nvec
					penTVec   = tvec
					penElLen = 0.5*(gg+1)*lenElem[im] + \
					           0.5*(1-gg)*lenElem[im-1]
					

					#print 'normVec - = ', normVec[penID]
					#print 'penVec    = ', penVec
					#print 'normVec + = ', normVec[penID+1]
					#print 'penDepth  = ', penDepth
					#print 'penXi     = ', penXi
					#print '\n'
					
			elif lms <0.2*lenElem[im] and xi>=1.-alpha and xi<1.+epsXi and not im==(self.numMasterNodes-2):
				# Hermite interpolation with next element
				
				gg = sco.fsolve(self.getHermiteG, (xi-1.)/alpha, args=(alpha,sLoc,mLoc[im],mLoc[im+1],mLoc[im+2]) )[0]
				dep,nvec,tvec = self.getHermiteDepth(gg,alpha,sLoc,mLoc[im],mLoc[im+1],mLoc[im+2])

				if gg<-1 or gg>1:
					print 'WARNING: gg = ', gg
				
				if dep < 0.:
					contactFound = True
					penDepth = -dep
					penID    = im           if gg<0 else im+1
					penXi    = alpha*gg+1.  if gg<0 else alpha*gg
					penNVec   = nvec
					penTVec   = tvec
					penElLen = 0.5*(gg+1)*lenElem[im+1] + \
					           0.5*(1-gg)*lenElem[im]

					#print 'normVec - = ', normVec[penID]
					#print 'penVec    = ', penVec
					#print 'normVec + = ', normVec[penID+1]
					#print 'penDepth  = ', penDepth
					#print 'penXi     = ', penXi
					#print '\n'
					
			# Increment im
			im+=1

		return (contactFound, {'penDepth':penDepth,'penID':penID,'penXi':penXi,'penNVec':penNVec,'penTVec':penTVec,'penElLen':penElLen})
	"""
	cpdef getPenalty(self,double penDepth):
		#return (self.penaltyParameter*penDepth**3,3*self.penaltyParameter*penDepth**2)
		return (self.penaltyParameter*penDepth**2,2*self.penaltyParameter*penDepth)
		#return (self.penaltyParameter*penDepth,self.penaltyParameter)
