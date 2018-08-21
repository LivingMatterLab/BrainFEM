# -*- coding: utf-8 -*-
# In this mechanism, elements behave according to the dynein mechanism:
# contract, release, extend, reattach, contract, etc.
from Mechanism cimport *
from CBARX cimport *
from ModelContainer cimport *
from Solver cimport *

import numpy as np
cimport numpy as np
from math import *
cimport OutputHelper as oh

cdef class MECH_EL03(Mechanism):
	def __init__(self, tCont, tDest, tCrea, minInitStretch, minStretch, maxInitStretch, maxStretch,activeStretch):
		super(MECH_EL03,self).__init__();

		self.tCont 			= tCont
		self.tDest 			= tDest
		self.tCrea 			= tCrea
		self.minInitStretch = minInitStretch
		self.minStretch 	= minStretch
		self.maxInitStretch = maxInitStretch
		self.maxStretch 	= maxStretch

		self.activeStretch 	= activeStretch

	cpdef Initialize(self, Element el, ModelContainer mc):
		cdef bint MT_active
		cdef double randMT


		# First check whether this element has active MT on both sides
		MT_active = True;
		try:
			if ( (el.nodes[0].elMinus.state == MicrotubuleInactive and el.nodes[0].elPlus.state == MicrotubuleInactive) or \
				 (el.nodes[1].elMinus.state == MicrotubuleInactive and el.nodes[1].elPlus.state == MicrotubuleInactive)):
				MT_active = False;
		except:
			MT_active = False;
			
		# Check whether crossLink should be attached or detached (50-50 chance)
		# This is 50-50 chance, and based whether MT elements are active
		randMT = np.random.rand()
		if(randMT<0.5*self.tCont/(self.tCont+self.tDest) and MT_active):
			# Attached slack
			el.dummyNodes = []
			el.timeToNextEvent = np.random.rand()*self.tCont
			el.state           = MotorAttachedSlack
		elif(randMT<0.5 and MT_active):
			# Attached taut
			el.dummyNodes = []
			el.restLength*=self.activeStretch
			el.timeToNextEvent = np.random.rand()*self.tDest
			el.state           = MotorAttachedTaut
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
		cdef bint is_Rel, is_Active, is_to_be_Rel, firstTry
		cdef np.ndarray[np.int_t, ndim=1] dof

		nCreate = 0

		if el.state==MotorJustAttachedSlack:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorAttachedSlack;
		if el.state==MotorJustAttachedTaut:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorAttachedTaut;
		if el.state==MotorJustReleased:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorReleased;

		if(el.timeToNextEvent<=0):
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
			
			if el.state==MotorAttachedSlack or el.state==MotorJustAttachedSlack:
				# Contract element by changing its restlength
				el.restLength*=self.activeStretch
				el.timeToNextEvent = np.random.rand()*self.tDest
				el.state           = MotorJustAttachedTaut
			elif el.state==MotorAttachedTaut or el.state==MotorJustAttachedTaut:
				# Destroy element by changing its connectivity to the storage nodes of the model
				el.dummyNodes = el.nodes
				el.nodes = mc.storageNodes
				el.restLength/=self.activeStretch
				el.timeToNextEvent = np.random.rand()*self.tCrea
				el.state           = MotorJustReleased
				nCreate -=1
			elif el.state==MotorReleased or el.state==MotorJustReleased:
				# Create element by changing its connectivity to its dummyNodes
				el.nodes = el.dummyNodes
				el.dummyNodes = []
				el.timeToNextEvent = np.random.rand()*self.tCont
				el.state           = MotorJustAttachedSlack
				nCreate +=1

				# First check whether the node has active MT on both sides
				MT_active = True;
				try:
					if ( (el.nodes[0].elMinus.state == MicrotubuleInactive and el.nodes[0].elPlus.state == MicrotubuleInactive) or \
						 (el.nodes[1].elMinus.state == MicrotubuleInactive and el.nodes[1].elPlus.state == MicrotubuleInactive)):
						MT_active = False;
				except:
					MT_active = False;

				if MT_active:
					# Change one of the element nodes, such that element length is closest 
					# to rest length of the element
					eidMove = 1
						
					# We can move, first compute initial distance
					dist0 = el.CurrentLength(s.dc)

					# Dynein always moves in minus direction 
					moveDir = -1;
					dist  = el.CurrentLength(s.dc)

					# Move in the computed direction until we are connected to the node that 
					# gives an element length that is closest to the rest length.
					while dist<el.restLength and not el.nodes[eidMove].elMinus==None:
						dist0 = dist;
						el.nodes[eidMove] = el.nodes[eidMove].elMinus.nodes[0]
						dist  = el.CurrentLength(s.dc)
						frac = (el.restLength-dist0)/(dist-dist0)

					# Check whether we moved 1 element too far
					if frac<0.5 and frac>-0.5 and not el.nodes[eidMove].elPlus==None:
						dist0 = dist;
						# We moved to far, so undo last move
						el.nodes[eidMove] = el.nodes[eidMove].elPlus.nodes[1]
						dist  = el.CurrentLength(s.dc)
						frac = (el.restLength-dist0)/(dist-dist0)

					
					if(dist/el.restLength>self.maxInitStretch or dist/el.restLength<self.minInitStretch ):
						# If not, recreate this cross-link randomly along axon
						i = 0;
						iMax = 100;
						while(i<iMax and (dist/el.restLength>self.maxInitStretch or dist/el.restLength<self.minInitStretch)):
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
						"""
						# If not, recreate this cross-link randomly along axon
						i = 0;
						iMax = 100;
						while(i<iMax and (dist/el.restLength>self.maxInitStretch or dist/el.restLength<self.minInitStretch)):
							i+=1
							
							# New x location of nodes (in reference coordinates)
							nid = i%2
							if i<5:
								xDest = el.nodes[nid].x
							else:
								xDest = np.random.rand()*mc.lAxon
								
								# Move first node
								if xDest<el.nodes[nid].x:	
									# Move to minus side
									while(xDest<el.nodes[nid].x and not el.nodes[nid].elMinus==None):
										el.nodes[nid] = el.nodes[nid].elMinus.nodes[0]
								else:
									# Move to plus side
									while(xDest>el.nodes[nid].x and not el.nodes[nid].elPlus==None):
										el.nodes[nid] = el.nodes[nid].elPlus.nodes[1]

							# Move second node
							nid = (nid+1)%2
							dist  = el.CurrentLength(s.dc)
							if xDest<el.nodes[nid].x:	
								# Move to minus side
								while(dist/el.restLength>self.maxInitStretch and not el.nodes[nid].elMinus==None):
									el.nodes[nid] = el.nodes[nid].elMinus.nodes[0]
									dist  = el.CurrentLength(s.dc)
							if xDest>el.nodes[nid].x:	
								# Move to plus side
								while(dist/el.restLength>self.maxInitStretch and not el.nodes[nid].elPlus==None):
									el.nodes[nid] = el.nodes[nid].elPlus.nodes[1]
									dist  = el.CurrentLength(s.dc)
						"""
						# Check which end has to be walking domain
						frac0 = 1.; frac1=1.
						for mtlID in range(len(mc.MTs)):
							for mtID in range(len(mc.MTs[mtlID])):
								n0 = mc.MTs[mtlID][mtID].n0
								n1 = mc.MTs[mtlID][mtID].n1
								n2 = mc.MTs[mtlID][mtID].n2
								if(el.nodes[0].localID>=n0 and el.nodes[0].localID<=n2):
									frac0 = 1.*(el.nodes[0].localID-n0)/(n1-n0)
								if(el.nodes[1].localID>=n0 and el.nodes[1].localID<=n2):
									frac1 = 1.*(el.nodes[1].localID-n0)/(n1-n0)

						# For dynein motors, the walking side is towards the minus end of the MT
						# The walking side is the second node the element. Switch order of nodes 
						# if necessary
						if (mc.optionDynein==0 and frac0<frac1) or (mc.optionDynein==2 and frac0>frac1):
							el.nodes = [el.nodes[1],el.nodes[0]]
						
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
		if((el.state==MotorAttachedSlack or el.state==MotorJustAttachedSlack  or \
			el.state==MotorAttachedTaut or el.state==MotorJustAttachedTaut)\
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
			print "Destroyed: MaxStretch"

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
			
			
			"""
			print "Destroyed: Inactive MT"
			print el.nodes[0]
			print el.nodes[1]
			try:
				print el.nodes[0].elMinus.state, el.nodes[1].elMinus.state, el.nodes[0].elPlus.state, el.nodes[1].elPlus.state
			except:
				pass
			"""

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