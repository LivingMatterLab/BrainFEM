# -*- coding: utf-8 -*-
# In this mechanism, elements behave according to the myosin mechanism:
# one (randomly picked) end moves towards plus end (no detachments)
from Mechanism cimport *
from CBARX cimport *
from ModelContainer cimport *
from Solver cimport *

import numpy as np
cimport numpy as np
from math import *
cimport OutputHelper as oh

cdef class MECH_EL05(Mechanism):
	def __init__(self, tMove, minStretch, maxStretch,activeStretch):
		super(MECH_EL05,self).__init__();

		self.tMove 			= tMove
		self.minStretch 	= minStretch
		self.maxStretch 	= maxStretch

		self.activeStretch 	= activeStretch

	cpdef Initialize(self, Element el, ModelContainer mc):
		cdef bint MT_active
		cdef double randMT
		# Check whether crossLink should be attached or detached
		# This is based on whether MT elements are active
		MT_active = True;
		try:
			if ( (el.nodes[0].elMinus.state == MicrotubuleInactive and el.nodes[0].elPlus.state == MicrotubuleInactive) or \
				 (el.nodes[1].elMinus.state == MicrotubuleInactive and el.nodes[1].elPlus.state == MicrotubuleInactive)):
				MT_active = False;
		except:
			MT_active = False;
			
		if MT_active:
			# Attached slack
			el.dummyNodes = []
			el.timeToNextEvent = (np.random.rand()-0.5)*self.tMove # -0.5 to initialize half of crosslinks to be stretched
			el.state           = MotorAttachedSlack
		else:
			# Detached
			el.dummyNodes = el.nodes
			el.nodes = mc.storageNodes
			el.timeToNextEvent = 1.e9
			el.state           = MotorReleased

	cpdef Apply(self, Element el, ModelContainer mc, Solver s):
		cdef int i
		cdef int nCreate, ndof, countK, eidMove
		cdef double dist, dist0, frac, xDest
		cdef bint is_Rel, is_Active, is_to_be_Rel, firstTry
		cdef np.ndarray[np.int_t, ndim=1] dof

		nCreate = 0

		if el.state==MotorJustAttachedTaut:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorAttachedSlack;
		if el.state==MotorJustReleased:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorReleased;

		if(el.timeToNextEvent<=0):
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
			
			# Myosin will always attempt to stretch itself to 
			# its activeStretch by walking to the plus end of 
			# the actin filament
			if el.state==MotorAttachedSlack or el.state==MotorJustAttachedTaut:
				# Create element by changing its connectivity to its dummyNodes
				el.timeToNextEvent = np.random.rand()*self.tMove
				el.state           = MotorJustAttachedTaut

				# Change one of the element nodes, such that element length is closest 
				# to rest length of the element
				eidMove = int(2*np.random.rand())
					
				# We can move, first compute initial distance
				dist0 = el.CurrentLength(s.dc)

				# Myosin always moves in plus direction 
				moveDir = 1;
				dist  = el.CurrentLength(s.dc)

				# Move in plus direction until we are connected to the node that 
				# gives an element length that is closest to the rest length.
				i = 0; iMax = 100;
				while i<iMax and (dist<el.restLength*self.activeStretch or dist<=dist0) and not el.nodes[eidMove].elPlus==None:
					i+=1;
					dist0 = dist;
					plusID = el.nodes[eidMove].elPlus.plusID;
					el.nodes[eidMove] = el.nodes[eidMove].elPlus.nodes[plusID]
					dist  = el.CurrentLength(s.dc)
					frac = 0.5 if np.abs(dist-dist0)<1.e-9 else (el.restLength*self.activeStretch-dist0)/(dist-dist0)

				# Check whether we moved 1 element too far
				if frac<0.5 and frac>-0.5 and not el.nodes[eidMove].elMinus==None:
					dist0 = dist;
					# We moved to far, so undo last move
					minusID = el.nodes[eidMove].elMinus.minusID;
					el.nodes[eidMove] = el.nodes[eidMove].elMinus.nodes[minusID]
					dist  = el.CurrentLength(s.dc)
					frac = (el.restLength*self.activeStretch-dist0)/(dist-dist0)


				if i==iMax:
					# If not converged, recreate this cross-link randomly along axon
					i = 0;
					iMax = 100;
					while(i<iMax and (dist/el.restLength>self.maxStretch or dist/el.restLength<self.minStretch)):
						i+=1
						
						# New x location of nodes (in reference coordinates)
						xDest = np.random.rand()*mc.lAxon

						# Move first node
						nid = i%2
						if xDest<el.nodes[nid].x:
							# Move to negative x direction	
							while(xDest<el.nodes[nid].x and not \
								  (el.nodes[nid].elMinus==None or el.nodes[nid].elMinus==None)):
								el.nodes[nid] = mc.nodes[el.nodes[nid].localID-1];
						else:
							# Move to positive x direction	
							while(xDest>el.nodes[nid].x and not \
								  (el.nodes[nid].elMinus==None or el.nodes[nid].elMinus==None)):
								el.nodes[nid] = mc.nodes[el.nodes[nid].localID+1];

						# Move second node
						nid = (nid+1)%2
						dist  = el.CurrentLength(s.dc)
						if xDest<el.nodes[nid].x:	
							# Move to negative x direction	
							while(dist/el.restLength>self.maxStretch and not \
								  (el.nodes[nid].elMinus==None or el.nodes[nid].elMinus==None)):
								el.nodes[nid] = mc.nodes[el.nodes[nid].localID-1]
								dist  = el.CurrentLength(s.dc)
						if xDest>el.nodes[nid].x:	
							# Move to positive x direction	
							while(dist/el.restLength>self.maxStretch and not \
								  (el.nodes[nid].elMinus==None or el.nodes[nid].elMinus==None)):
								el.nodes[nid] = mc.nodes[el.nodes[nid].localID+1]
								dist  = el.CurrentLength(s.dc)

					
					# If attempts to recreate failed, then destroy element again
					if i==iMax:
						print "Destroyed: InitStretch: ", dist/el.restLength
						el.dummyNodes = el.nodes
						el.nodes = mc.storageNodes
						el.timeToNextEvent = 1.e9
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
		if((el.state==MotorAttachedSlack or el.state==MotorJustAttachedSlack  or \
			el.state==MotorAttachedTaut or el.state==MotorJustAttachedTaut)\
			and el.CurrentLength(s.dc)/el.restLength>self.maxStretch):
			# Add element to storage
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
			
			# Destroy element by changing its connectivity to the storage nodes of the model
			el.dummyNodes = el.nodes
			el.nodes = mc.storageNodes
			el.timeToNextEvent = 1.e9
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
			el.timeToNextEvent = 1.e9
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