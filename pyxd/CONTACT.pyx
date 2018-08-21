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
		normVec = []
		for im in range(self.numMasterNodes-1):
			ex = mLoc[im+1]-mLoc[im]
			ex = ex/np.linalg.norm(ex)
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
			for im in range(self.numMasterNodes-1):
				# Compute vector from center element to slave node
				ems  = sLoc-0.5*(mLoc[im] +mLoc[im+1])
				lms  = np.inner(normVec[im],ems)

				# Compute xi, indicating whether slave node is 'inside' two master nodes
				xi = np.inner(sLoc-mLoc[im],mLoc[im+1]-mLoc[im])/(np.linalg.norm(mLoc[im+1]-mLoc[im])**2)

				if lms <0 and np.abs(lms)<penDepth and xi>0. and xi<1.:
					penDepth = -lms
					penID = im

			# Check if penetration occurred and act accordingly
			if penDepth < np.inf:
				vid_M = (self.numSlaveNodes+penID)*ndim
				penForce =self.penaltyParameter*penDepth**2

				kid = [i for i in range(vid_S,vid_S+ndim)]
				kid.extend([i for i in range(vid_M,vid_M+2*ndim)])


				nn = np.outer(normVec[penID],normVec[penID])
				mat[np.ix_(kid,kid)]+= 2*self.penaltyParameter*penDepth*np.bmat([[ nn, -.5* nn, -.5* nn],[-.5* nn, .5* nn,np.zeros((ndim,ndim))],\
					                                      [-.5* nn, np.zeros((ndim,ndim)),.5* nn]])
			
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
		normVec = []
		for im in range(self.numMasterNodes-1):
			ex = mLoc[im+1]-mLoc[im]
			ex = ex/np.linalg.norm(ex)
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
			for im in range(self.numMasterNodes-1):
				# Compute vector from center element to slave node
				ems  = sLoc-0.5*(mLoc[im] +mLoc[im+1])
				lms  = np.inner(normVec[im],ems)

				# Compute xi, indicating whether slave node is 'inside' two master nodes
				xi = np.inner(sLoc-mLoc[im],mLoc[im+1]-mLoc[im])/(np.linalg.norm(mLoc[im+1]-mLoc[im])**2)
				
				if lms <0 and np.abs(lms)<penDepth and xi>0. and xi<1.:
					penDepth = -lms
					penID = im


			# Check if penetration occurred and act accordingly
			if penDepth < np.inf:
				vid_M = (self.numSlaveNodes+penID)*ndim
				penForce = self.penaltyParameter*penDepth**2

				# Add penForce to loadVec
				vec[vid_S:vid_S+ndim]+=  penForce*normVec[penID]
				
				vec[vid_M:vid_M+ndim]        -=  .5*penForce*normVec[penID]
				vec[vid_M+ndim:vid_M+2*ndim] -=  .5*penForce*normVec[penID]
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
