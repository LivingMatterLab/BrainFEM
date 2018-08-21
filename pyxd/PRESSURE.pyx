# -*- coding: utf-8 -*-
import numpy as np
cimport numpy as np
import scipy.optimize as sco

cdef class PRESSURE(object):
	def __init__(self, nodes, pressure, amplitude):
		self.nodes       = nodes
		self.numNodes    = len(nodes)
		self.pressure    = pressure
		self.amplitude   = amplitude

	def __str__(self):
		return "PRESSURE: \t nodes = ["+ ', '.join(str(x.localID) for x in self.nodes)+"]" \
		               +"\t pressure = "+ str(self.pressure)


	cpdef BuildMatrix(self, DataContainer dc):
		cdef int ndim, vidN, nID
		cdef double press, lenElem
		cdef np.ndarray[double,ndim=1] ex,normVec
		cdef np.ndarray[double,ndim=2] nn,Inn,mat
		cdef list nLoc
		

		# Current pressure
		press = self.pressure*self.amplitude.Get(dc.time)

		ndim = 2;
		mat = np.zeros((self.datRID[1]-self.datRID[0],self.datRID[1]-self.datRID[0]));
		
		# Collect current locations of all nodes
		nLoc = [];
		for n in self.nodes:
			nLoc.append(np.array(n.loc)+np.array(n.DispDof(dc)))
		
		# Compute current normal vectors to each master line element
		vidN = -ndim;
		for eID in range(self.numNodes-1):
			# normVec and lenElem
			ex       = nLoc[eID+1]-nLoc[eID]
			lenElem  = np.linalg.norm(ex)
			ex       = ex/lenElem
			normVec  = np.array([-ex[1],ex[0]])
			nn = np.outer(normVec,normVec)
			nt = np.outer(normVec,ex)
			Inn = np.eye(len(normVec))-nn
			
			# nodal id's of this element
			vidN +=ndim
			vid   = [i for i in range(vidN,vidN+2*ndim)]

			# Add to mat
			"""
			print 'nn = \n', nn
			print 'nt = \n', nt
			print 'Inn = \n', Inn
			print '\n'
			"""
			mat[np.ix_(vid,vid)]-= 0.5*press*np.bmat([[-nt, nt],[-nt,nt]])
			#mat[np.ix_(vid,vid)]-= 0.5*press*np.bmat([[-nt.T, -nt.T],[nt.T,nt.T]])
			
			mat[np.ix_(vid,vid)]-= 0.5*press*np.bmat([[-Inn,Inn],[-Inn,Inn]])/lenElem
			

			#nn = np.array([[0.,-1.],[1.,0.]])
			#mat[np.ix_(vid,vid)]-= 0.5*press*np.bmat([[-nn, -nn],[nn,nn]]).T

		return mat

	cpdef BuildVector(self, DataContainer dc):
		cdef int ndim, vidN, nID
		cdef double press, lenElem
		cdef np.ndarray[double,ndim=1] vec, ex,normVec
		cdef list nLoc

		# Current pressure
		press = self.pressure*self.amplitude.Get(dc.time)

		ndim = 2;
		vec = np.zeros(self.datRID[1]-self.datRID[0]);
		
		# Collect current locations of all nodes
		nLoc = [];
		for n in self.nodes:
			nLoc.append(np.array(n.loc)+np.array(n.DispDof(dc)))
		
		# Compute current normal vectors to each master line element
		vidN = -ndim;
		for eID in range(self.numNodes-1):
			# normVec and lenElem
			ex       = nLoc[eID+1]-nLoc[eID]
			lenElem  = np.linalg.norm(ex)
			ex       = ex/lenElem
			normVec  = np.array([-ex[1],ex[0]])
			
			# nodal id's of this element
			vidN +=ndim
			vid   = [i for i in range(vidN,vidN+2*ndim)]

			# Add to vec
			#print lenElem
			#print press
			#print 0.5*press*lenElem*np.concatenate([-normVec,-normVec])
			vec[vid]-= 0.5*press*lenElem*np.concatenate([normVec,normVec])

		return vec
