# -*- coding: utf-8 -*-
import numpy as np
cimport numpy as np

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
			vid_S+=ndim
			sLoc  = np.array(s.loc)+np.array(s.DispDof(dc))
			sLoc0 = np.array(s.loc)+np.array(s.DispDof0(dc))
			penDepth = np.inf
			penID    = -1
			penXi    = 0.5
			for im in range(self.numMasterNodes-1):
				# Compute xi, indicating whether slave node is 'inside' two master nodes
				xi = np.inner(sLoc-mLoc[im],tangVec[im])/lenElem[im]

				# Compute vector from center element to slave node
				#ems  = sLoc-0.5*(mLoc[im] +mLoc[im+1])
				ems  = sLoc-(1-xi)*mLoc[im]-xi*mLoc[im+1]
				lms  = np.inner(normVec[im],ems)
				
				if lms <0 and np.abs(lms)<penDepth and xi>0. and xi<1.:
					penDepth = -lms
					penID = im
					penXi = xi

			# Check if penetration occurred and act accordingly
			if penDepth < np.inf:
				vid_M = (self.numSlaveNodes+penID)*ndim
				penForce =self.penaltyParameter*penDepth**2


				# Compute normal/tangent vectors
				Ns  = np.concatenate([normVec[penID],-(1.-penXi)*normVec[penID],-penXi*normVec[penID]])
				N0s = np.concatenate([np.zeros(ndim),           -normVec[penID],       normVec[penID]])
				Ts  = np.concatenate([tangVec[penID],-(1.-penXi)*tangVec[penID],-penXi*tangVec[penID]])
				T0s = np.concatenate([np.zeros(ndim),           -tangVec[penID],       tangVec[penID]])
				
				# Add penForce to vec
				kid = [i for i in range(vid_S,vid_S+ndim)]
				kid.extend([i for i in range(vid_M,vid_M+2*ndim)])

				Kd = -(0*np.outer(N0s,Ts)+0*np.outer(Ts,N0s) +\
					   penDepth/lenElem[penID]*np.outer(N0s,N0s))/lenElem[penID]

				"""
				print 'Kd=\n',penForce*Kd
				print 'NN=\n',2*self.penaltyParameter*penDepth*np.outer(Ns,Ns)
				print '\n'
				"""

				mat[np.ix_(kid,kid)]+= (2*self.penaltyParameter*penDepth*np.outer(Ns,Ns)+penForce*Kd)
			
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
			vid_S+=ndim
			sLoc  = np.array(s.loc)+np.array(s.DispDof(dc))
			sLoc0 = np.array(s.loc)+np.array(s.DispDof0(dc))
			penDepth = np.inf
			penID    = -1
			penXi    = 0.5
			for im in range(self.numMasterNodes-1):
				# Compute xi, indicating whether slave node is 'inside' two master nodes
				xi = np.inner(sLoc-mLoc[im],tangVec[im])/lenElem[im]

				# Compute vector from center element to slave node
				#ems  = sLoc-0.5*(mLoc[im] +mLoc[im+1])
				ems  = sLoc-(1-xi)*mLoc[im]-xi*mLoc[im+1]
				lms  = np.inner(normVec[im],ems)

				if lms <0 and np.abs(lms)<penDepth and xi>0. and xi<1.:
					penDepth = -lms
					penID = im
					penXi = xi


			# Check if penetration occurred and act accordingly
			if penDepth < np.inf:
				#print 'Contact, penXi = ', penXi

				vid_M = (self.numSlaveNodes+penID)*ndim
				penForce = self.penaltyParameter*penDepth**2

				# Compute normal/tangent vectors
				Ns  = np.concatenate([normVec[penID],-(1.-penXi)*normVec[penID],-penXi*normVec[penID]])
				
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
