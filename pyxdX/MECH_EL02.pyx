# -*- coding: utf-8 -*-
# In this mechanism, elements reattach to the nodes that is results in
# an element length closest to its rest length
# It is then check that the element does not connect to an inactive mt
from Mechanism cimport *
from CBARX cimport *
from ModelContainer cimport *
from Solver cimport *

import numpy as np
cimport numpy as np
from math import *
cimport OutputHelper as oh

cdef class MECH_EL02(Mechanism):
	def __init__(self,tCrea, tDest, maxInitStretch, maxStretch):
		super(MECH_EL02,self).__init__();

		self.tCrea 			= tCrea
		self.tDest 			= tDest
		self.maxInitStretch = maxInitStretch
		self.maxStretch 	= maxStretch

	cpdef Initialize(self, Element el, ModelContainer mc):
		cdef bint MT_active
		# Check whether crossLink should be attached or detached (50-50 chance)
		# This is 50-50 chance, and based whether MT elements are active
		MT_active = True;
		try:
			if ( (el.nodes[0].elMinus.state == MicrotubuleInactive and el.nodes[0].elPlus.state == MicrotubuleInactive) or \
				 (el.nodes[1].elMinus.state == MicrotubuleInactive and el.nodes[1].elPlus.state == MicrotubuleInactive)):
				MT_active = False;
		except:
			MT_active = False;
			
		if(np.random.rand()<0.5 and MT_active):
			# Attached
			el.dummyNodes = []
			el.timeToNextEvent = np.random.rand()*self.tDest
			el.state           = MotorAttachedSlack
		else:
			# Detached
			el.dummyNodes = el.nodes
			el.nodes = mc.storageNodes
			el.timeToNextEvent = np.random.rand()*self.tCrea
			el.state           = MotorReleased

	cpdef Apply(self, Element el, ModelContainer mc, Solver s):
		cdef int i
		cdef int nCreate, ndof, countK, eidMove
		cdef double dist, dist0, frac, xDest
		cdef bint is_Rel, is_Active, is_to_be_Rel
		cdef np.ndarray[np.int_t, ndim=1] dof

		nCreate = 0

		if el.state==MotorJustAttachedSlack:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorAttachedSlack;
		if el.state==MotorJustReleased:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorReleased;

		if(el.timeToNextEvent<=0):
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
			
			if el.state==MotorAttachedSlack or el.state==MotorJustAttachedSlack:
				# Destroy element by changing its connectivity to the storage nodes of the model
				el.dummyNodes = el.nodes
				el.nodes = mc.storageNodes
				el.timeToNextEvent = np.random.rand()*self.tCrea
				el.state           = MotorJustReleased
				nCreate -=1
			elif el.state==MotorReleased or el.state==MotorJustReleased:
				# Create element by changing its connectivity to its dummyNodes
				el.nodes = el.dummyNodes
				el.dummyNodes = []
				el.timeToNextEvent = np.random.rand()*self.tDest
				el.state           = MotorJustAttachedSlack
				nCreate +=1

				# Change one of the element nodes, such that element length is closest 
				# to rest length of the element
				eidMove = np.random.randint(0,2)

				"""
				# First check whether our moving node has a plus and minus element
				# If not, recreate this cross-link randomly along axon
				i = 0;
				while( (el.nodes[eidMove].elPlus==None or el.nodes[eidMove].elMinus==None) and i<10):
					i+=1
					if el.nodes[eidMove].elPlus==None:
						# New x location of node (in reference coordinates)
						xDest = np.random.rand()*el.nodes[eidMove].x
						# Move
						while(xDest<el.nodes[eidMove].x and not el.nodes[eidMove].elMinus.nodes[0].elMinus==None):
							el.nodes[eidMove] = el.nodes[eidMove].elMinus.nodes[0]
					else:
						# New x location of node (in reference coordinates)
						xDest = el.nodes[eidMove].x + np.random.rand()*(mc.lAxon - el.nodes[eidMove].x)
						# Move
						while(xDest>el.nodes[eidMove].x and not el.nodes[eidMove].elPlus.nodes[1].elPlus==None):
							el.nodes[eidMove] = el.nodes[eidMove].elPlus.nodes[1]

					# Set other node to move
					eidMove = (eidMove+1)%2
				"""
				# We can move, first compute initial distance
				dist0 = el.CurrentLength(s.dc)

				# Guess we have to move in positive x direction
				try:
					el.nodes[eidMove] = mc.nodes[el.nodes[eidMove].localID+1];
					dist  = el.CurrentLength(s.dc)
					frac = 0.5 if np.abs(dist-dist0)<1.e-9 else (el.restLength-dist0)/(dist-dist0)		
					if frac>0:
						# Our guess was right, we move in positive x direction
						moveDir = 1;
					else:
						raise 0  # Will go to except statement
				except:
					# Our guess was wrong, we have to move in negative x direction
					moveDir = -1;
					try:
						el.nodes[eidMove] = mc.nodes[el.nodes[eidMove].localID-1];
						dist  = el.CurrentLength(s.dc)
						frac = 0.5 if np.abs(dist-dist0)<1.e-9 else (el.restLength-dist0)/(dist-dist0)
					except:
						frac = 0.5

				# Move in the computed direction until we are connected to the node that 
				# gives an element length that is closest to the rest length.
				i = 0; iMax = 100;
				while i<iMax and (dist<el.restLength or dist<=dist0) and not \
								  (el.nodes[eidMove].elMinus==None or el.nodes[eidMove].elPlus==None):
					i+=1;
					dist0 = dist;
					el.nodes[eidMove] = mc.nodes[el.nodes[eidMove].localID+moveDir];
					dist  = el.CurrentLength(s.dc)
					frac = 0.5 if np.abs(dist-dist0)<1.e-9 else (el.restLength-dist0)/(dist-dist0)

				# Check whether we moved 1 element too far
				if frac<0.5 and frac>-0.5 and not \
					(el.nodes[eidMove].elMinus==None or el.nodes[eidMove].elPlus==None):
					el.nodes[eidMove] = mc.nodes[el.nodes[eidMove].localID-moveDir];
					dist  = el.CurrentLength(s.dc)
				
				# Check crosslink is attached correctly, otherwise, recreate crosslink randomly along axon
				if(i==iMax or dist/el.restLength>self.maxInitStretch):
					# If not, recreate this cross-link randomly along axon
					i = 0;
					iMax = 100;
					while(i<iMax and (dist/el.restLength>self.maxInitStretch)):
						i+=1

						# Move first node
						nid = i%2
						if i<5:
							xDest = el.nodes[nid].x
						else:
							xDest = np.random.rand()*mc.lAxon
							
							if xDest<el.nodes[nid].x:
								# Move to negative x direction	
								while(xDest<el.nodes[nid].x and not \
									  (el.nodes[nid].elMinus==None or el.nodes[nid].elPlus==None)):
									el.nodes[nid] = mc.nodes[el.nodes[nid].localID-1];
							else:
								# Move to positive x direction	
								while(xDest>el.nodes[nid].x and not \
									  (el.nodes[nid].elMinus==None or el.nodes[nid].elPlus==None)):
									el.nodes[nid] = mc.nodes[el.nodes[nid].localID+1];
							
						# Move second node
						nid = (nid+1)%2
						dist  = el.CurrentLength(s.dc)
						if xDest<el.nodes[nid].x:	
							# Move to negative x direction	
							while(dist/el.restLength>self.maxInitStretch and not \
								  (el.nodes[nid].elMinus==None or el.nodes[nid].elPlus==None)):
								el.nodes[nid] = mc.nodes[el.nodes[nid].localID-1]
								dist  = el.CurrentLength(s.dc)
						if xDest>el.nodes[nid].x:	
							# Move to positive x direction	
							while(dist/el.restLength>self.maxInitStretch and not \
								  (el.nodes[nid].elMinus==None or el.nodes[nid].elPlus==None)):
								el.nodes[nid] = mc.nodes[el.nodes[nid].localID+1]
								dist  = el.CurrentLength(s.dc)
				
					# If attempts to recreate failed, then destroy element again
					if i==iMax:
						print "Destroyed: InitStretch: ", dist/el.restLength
						el.dummyNodes = el.nodes
						el.nodes = mc.storageNodes
						el.timeToNextEvent = np.random.rand()*self.tCrea
						el.state           = MotorJustReleased
						nCreate -=1

			else:
				print "ERROR: el.state ", el.state, " does not have a next event."

			# Update rowR, rowK, colK in the solver
			dof = el.DofID();
			ndof = len(dof);

			countK = el.datKID[0]
			for i in range(ndof):
				s.rowR[el.datRID[0]+i]=dof[i]
				for j in range(ndof):
					s.rowK[countK] = dof[j]
					s.colK[countK] = dof[i]
					countK +=1

		# Check if cross link stretch is too large
		if((el.state==MotorAttachedSlack or el.state==MotorJustAttachedSlack)\
			and el.CurrentLength(s.dc)/el.restLength>self.maxStretch):
			# Add element to storage
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
			
			# Destroy element by changing its connectivity to the storage nodes of the model
			el.dummyNodes = el.nodes
			el.nodes = mc.storageNodes
			el.timeToNextEvent = np.random.rand()*self.tCrea
			el.state           = MotorJustReleased
			nCreate -=1

			# Update rowR, rowK, colK in the solver
			dof = el.DofID();
			ndof = len(dof);

			countK = el.datKID[0]
			for i in range(ndof):
				s.rowR[el.datRID[0]+i]=dof[i]
				for j in range(ndof):
					s.rowK[countK] = dof[j]
					s.colK[countK] = dof[i]
					countK +=1

		# Check if element is cross link attached to inactive MT
		is_Rel = False			# Element is a released cross-link
		is_Active = False;		# Element is connected to active MT
		is_to_be_Rel = False;	# Element needs to be released

		# Check if element is a released cross link
		if(el.state==MotorReleased or el.state==MotorJustReleased):
			is_Rel = True;
		# Active cross link if:
		# -  element is not released
		# -  element nodes have no empty elMinus or elPlus
		# -  element nodes have no inactive MT on both sides
		if( not is_Rel and not (\
			el.nodes[0].elMinus==None or el.nodes[0].elPlus==None or \
		   	el.nodes[1].elMinus==None or el.nodes[1].elPlus==None) and not (\
			(el.nodes[0].elMinus.state == MicrotubuleInactive and el.nodes[0].elPlus.state == MicrotubuleInactive) or \
			(el.nodes[1].elMinus.state == MicrotubuleInactive and el.nodes[1].elPlus.state == MicrotubuleInactive))):
			is_Active = True;
		# Check if cross link needs to be released
		if(not is_Rel and not is_Active):
			is_to_be_Rel = True
		
		if is_to_be_Rel:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
			
			# Destroy element by changing its connectivity to the storage nodes of the model.
			el.dummyNodes = el.nodes
			el.nodes = mc.storageNodes
			el.timeToNextEvent = np.random.rand()*self.tCrea
			el.state           = MotorJustReleased
			nCreate -= 1

			# Update rowR, rowK, colK in the solver
			dof = el.DofID();
			ndof = len(dof);

			countK = el.datKID[0]
			for i in range(ndof):
				s.rowR[el.datRID[0]+i]=dof[i]
				for j in range(ndof):
					s.rowK[countK] = dof[j]
					s.colK[countK] = dof[i]
					countK +=1
		return nCreate