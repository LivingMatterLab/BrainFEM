# -*- coding: utf-8 -*-
# DAMAGE mechanism, elements will detach based on force and rate of force
# according to the Bell model. Reattachment is based on random distribution in
# time
from Mechanism cimport *
from CBARX cimport *
from ModelContainer cimport *
from Solver cimport *

import numpy as np
cimport numpy as np
from math import *
cimport OutputHelper as oh
cimport ElementHelper as eh

cdef class MECH_EL06(Mechanism):
	def __init__(self, tDest0, f_beta, maxInitStretch, maxStretch):
		super(MECH_EL06,self).__init__();

		self.tDest0 		= tDest0
		self.f_beta         = f_beta
		self.maxInitStretch = maxInitStretch
		self.maxStretch 	= maxStretch

		self.tCountDown     = 100000.		# Time to countdown after previous event

	cpdef Initialize(self, Element el, ModelContainer mc):
		cdef bint MT_active
		cdef double randMT
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
			el.timeToNextEvent = self.tCountDown # Detachment will be based on force (rate)
			el.state           = MotorAttachedSlack
		else:
			# Detached
			el.dummyNodes = el.nodes
			el.nodes = mc.storageNodes
			el.timeToNextEvent = float("inf")
			el.state           = MotorReleased

	cpdef Apply(self, Element el, ModelContainer mc, Solver s):
		cdef int i
		cdef int nCreate, ndof, countK, eidMove
		cdef str strWrite
		cdef double dist, dist0, frac, xDest,dtime
		cdef bint is_Rel, is_Active, is_to_be_Rel
		cdef np.ndarray[np.int_t, ndim=1] dof

		nCreate = 0

		if el.state==MotorJustAttachedSlack:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorAttachedSlack;
		if el.state==MotorJustReleased:
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			el.state = MotorReleased;

		# If element is attached, detachment will be based on force and force rate
		if el.state==MotorAttachedSlack or el.state==MotorJustAttachedSlack:

			# Compute current stretch and corresponding force
			stretch = el.CurrentLength(s.dc)/el.restLength
			M  = el.property.Piola1Stiffness(stretch,0,s.dc)
			force = M['P']

			# Compute average force rate
			dtime = self.tCountDown-el.timeToNextEvent
			if dtime-s.dc.dt<=0 or s.dc.time>mc.tLoad:
				force_rate = 0.;
			else:
				force_rate = force/(dtime-s.dc.dt)		# Average force rate

			if force_rate<1.e-6:
				# STATIC EVENT
				# Compute probability of detachment during time step
				prob_detach = (self.CdfStatBreakage(dtime,force)- \
							   self.CdfStatBreakage(dtime-s.dc.dt,force))/ \
							  (1.-self.CdfStatBreakage(dtime-s.dc.dt,force))
			else:
				# DYNAMIC EVENT 
				# Predict force at end of this step
				force1 = force+force_rate*s.dc.dt

				# Compute probability of detachment during time step
				prob_detach = (self.CdfDynBreakage(abs(force1),abs(force_rate))- \
							  self.CdfDynBreakage(abs(force) ,abs(force_rate)) )/\
							  (1.-self.CdfDynBreakage(abs(force) ,abs(force_rate)))
			# Check whether element will detach
			if np.random.rand()<prob_detach:
				mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
				
				# Print data of this detachment (for statistics only)
				strWrite= '{: <16d}'.format(el.localID) + \
		           '{: <16.6e}'.format(force) + \
		           '{: <16.6e}'.format(force_rate) + \
		           '{: <16.6f}'.format(dtime-s.dc.dt)
				oh.WriteToOutput(mc,'Bell.txt',strWrite)

				# Destroy element by changing its connectivity to the storage nodes of the model
				el.dummyNodes = el.nodes
				el.nodes = mc.storageNodes
				el.timeToNextEvent = float("inf")
				el.state           = MotorJustReleased
				nCreate -=1

				# Update rowR, rowK, colK in the solver
				eh.UpdateConnectivities(el,s)
		elif el.state==MotorReleased or el.state==MotorJustReleased:
			# Compute probability of detachment during time step
			prob_attach = self.CdfStatBreakage(s.dc.dt,0.)
			# Check whether element will detach
			if np.random.rand()<prob_attach:
				mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})

				# Create element by changing its connectivity to its dummyNodes
				el.nodes = el.dummyNodes
				el.dummyNodes = []
				el.timeToNextEvent = self.tCountDown	# Detachment will be based on force (rate)
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
					eidMove = np.random.randint(0,2)

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
						
					# We can move, first compute initial distance
					dist0 = el.CurrentLength(s.dc)

					# Guess we have to move in plus direction
					el.nodes[eidMove] = el.nodes[eidMove].elPlus.nodes[1]
					dist  = el.CurrentLength(s.dc)
					frac = 0.5 if np.abs(dist-dist0)<1.e-9 else (el.restLength-dist0)/(dist-dist0)		
					if frac>0:
						# Our guess was right, we move in positive direction
						moveDir = 1;
					else:
						# Our guess was wrong, we have to move in the minus direction
						moveDir = -1;
						el.nodes[eidMove] = el.nodes[eidMove].elMinus.nodes[0].elMinus.nodes[0]
						dist  = el.CurrentLength(s.dc)
						frac = 0.5 if np.abs(dist-dist0)<1.e-9 else (el.restLength-dist0)/(dist-dist0)
					# Move in the computed direction until we are connected to the node that 
					# gives an element length that is closest to the rest length.
					if moveDir == 1:
						while frac>1. and not el.nodes[eidMove].elPlus==None:
							dist0 = dist;
							el.nodes[eidMove] = el.nodes[eidMove].elPlus.nodes[1];
							dist  = el.CurrentLength(s.dc)
							frac = (el.restLength-dist0)/(dist-dist0)
					else:
						while frac>1. and not el.nodes[eidMove].elMinus==None:
							dist0 = dist;
							el.nodes[eidMove] = el.nodes[eidMove].elMinus.nodes[0]
							dist  = el.CurrentLength(s.dc)
							frac = (el.restLength-dist0)/(dist-dist0)

					# Check whether we moved 1 element too far
					if frac<0.5 and frac>-0.5:
						# We moved to far, so undo last move
						if moveDir == 1:
							el.nodes[eidMove] = el.nodes[eidMove].elMinus.nodes[0]
						else:
							el.nodes[eidMove] = el.nodes[eidMove].elPlus.nodes[1]

					if(el.CurrentLength(s.dc)/el.restLength>self.maxInitStretch):
						# Destroy element again
						el.dummyNodes = el.nodes
						el.nodes = mc.storageNodes
						el.timeToNextEvent = float("inf")
						el.state           = MotorJustReleased
						nCreate -=1
			
			# Update rowR, rowK, colK in the solver
			eh.UpdateConnectivities(el,s)
		else:
			print "ERROR: el.state ", el.state, " is not accounted for."


		# Check if cross link stretch is too large
		if((el.state==MotorAttachedSlack or el.state==MotorJustAttachedSlack)\
			and el.CurrentLength(s.dc)/el.restLength>self.maxStretch):
			# Add element to storage
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
			                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
			print "Destroyed: maxStretch = ", el.CurrentLength(s.dc)/el.restLength
			# Destroy element by changing its connectivity to the storage nodes of the model
			el.dummyNodes = el.nodes
			el.nodes = mc.storageNodes
			el.timeToNextEvent = float("inf")
			el.state           = MotorJustReleased
			nCreate -=1


			# Update rowR, rowK, colK in the solver
			eh.UpdateConnectivities(el,s)

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
			el.timeToNextEvent = float("inf")
			el.state           = MotorJustReleased
			nCreate -= 1

			# Update rowR, rowK, colK in the solver
			eh.UpdateConnectivities(el,s)
		return nCreate


	# Probability of crosslink breakage according to static Bell model
	cpdef CdfStatBreakage(self, double time, double force):
		cdef double k, k0
		k0 = 1./self.tDest0
		k  = k0*np.exp(force/self.f_beta)
		return 1.-np.exp(-time*k)

	# Probability of crosslink breakage according to dynamic Bell model
	cpdef PdfDynBreakage(self, double force, double force_rate):
		cdef double k, k0
		k0 = 1./self.tDest0
		k  = k0*np.exp(force/self.f_beta)
		return k/force_rate*np.exp(-self.f_beta/force_rate*(k-k0))

	cpdef CdfDynBreakage(self, double force, double force_rate):
		cdef double k, k0
		k0 = 1./self.tDest0
		k  = k0*np.exp(force/self.f_beta)
		return 1.-np.exp(-self.f_beta/force_rate*(k-k0))


