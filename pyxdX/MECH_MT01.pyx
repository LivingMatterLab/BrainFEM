# -*- coding: utf-8 -*-
# Mechanism that allows MT to polymerize, depolymerize and to be stationary
from Mechanism cimport *
from ModelContainer cimport *
from Node cimport *
from Element cimport *
from CBARX2 cimport *
from Solver cimport *

import numpy as np
cimport numpy as np
from math import *

cdef class MECH_MT01(Mechanism):
	def __init__(self,tMTpoly,tMTstat,tMTdepoly,MTpolyRate,MTdepolyRate,fracLMT=0.,goalNMT=9.5):
		super(MECH_MT01,self).__init__();

		self.tMTpoly 		= tMTpoly
		self.tMTstat 		= tMTstat
		self.tMTdepoly 		= tMTdepoly
		self.MTpolyRate     = MTpolyRate
		self.MTdepolyRate   = MTdepolyRate

		# fracLMT governs how much the total MT length is controlled (0-no control, 1-full control)
		self.fracLMT        = fracLMT		# Fraction determining control strength
		self.goalNMT        = goalNMT	    # Target number of mt per cross section
		print "self.goalNMT        =",self.goalNMT

	cpdef Initialize(self, MT mt, ModelContainer mc):
		cdef int eid
		cdef double randMT, xec, x1, t_cycle
		cdef Element el

		# Initialize MTstate and timeToNextMTstate
		randMT = np.random.rand()
		t_cycle = self.tMTpoly+self.tMTstat+self.tMTdepoly
		if(randMT< self.tMTpoly /t_cycle): # Polymerizing
			mt.state = Polymerizing
			mt.timeToNextEvent = np.random.rand()*self.tMTpoly
		elif(randMT<(self.tMTpoly+self.tMTstat)/t_cycle): # Stationary
			mt.state = Stationary
			mt.timeToNextEvent = np.random.rand()*self.tMTstat
		else: # Depolymerizing
			mt.state = Depolymerizing
			mt.timeToNextEvent = np.random.rand()*self.tMTdepoly
		if mt.mtPlus==None and np.abs(mc.nodes[mt.n1].x-mc.lAxon)<1.e-6:		
			# Last MT in line, does never polymerize or depolymerize
			mt.state = Stationary
			mt.timeToNextEvent = float("inf")
		if mt.mtMinus==None and mc.nodes[mt.n0].x<1.e-6:
			# First MT, if attached to wall, does never polymerize or depolymerize
			mt.state = Stationary
			mt.timeToNextEvent = float("inf")

		# Initialize state of elements in this mt
		x1 = mc.nodes[mt.n1].x
		for eid in range(mt.e0,mt.e2+1):
			el = mc.elements[eid]

			# First element in mt is always present
			if mt.e0==mt.e1:
				raise("ERROR: mt.e0=mt.e1")
			if eid==mt.e0:
				el.state = Microtubule
				el.timeToNextEvent = float("inf")
				continue

			# Compute initial state
			if(eid>mt.e1):
				el.state = MicrotubuleInactive
			else:
				el.state = Microtubule

			# Compute initTimeToNextEvent of element based on MTstate
			xec = 0.5*(el.nodes[0].x + el.nodes[1].x) # x at center of element
			if mt.state==Polymerizing:
				if el.state == MicrotubuleInactive:
					el.timeToNextEvent = (xec-x1)/self.MTpolyRate
				else:
					el.timeToNextEvent = float("inf")
			elif mt.state==Stationary:
				el.timeToNextEvent = float("inf")
			else:
				if el.state == MicrotubuleInactive:
					el.timeToNextEvent = float("inf")
				else:
					el.timeToNextEvent = (x1-xec)/self.MTdepolyRate
				

	cpdef Apply(self, MT mt, ModelContainer mc, Solver s):
		cdef Element el
		cdef int eid
		cdef double x1, xec, randMT

		if(mt.timeToNextEvent<=0):
			self.ChangeMTState(mt,mc,s)
		elif(mt.state==Stationary): # No event is happening in this MT
			pass;
		elif(mt.state==Polymerizing): #Event may happen in this MT
			mc.mtRestore.append({'mt':mt,'n1':mt.n1,'e1':mt.e1, 'state':mt.state, \
								    'timeToNextEvent':mt.timeToNextEvent})
			
			for eid in range(mt.e1,mt.e2+1): # Loop from - to + side 
				el = mc.elements[eid]
				if el.timeToNextEvent<=0:
					"""
					print "Before polymerizing..."
					print mt
					for nid in range(mt.n0,mt.n2+3):
						print mc.nodes[nid]
					for eid in range(mt.e0,mt.e2+1):
						print mc.elements[eid]
					"""
					mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})

					# Set new state and timeToNextEvent in this element
					el.state = Microtubule
					el.timeToNextEvent = float("inf")
					# Update n1 and e1 in this element
					mt.n1 +=1
					mt.e1 +=1

					# Update pointers
					el.nodes[0].elPlus = el 
					if not mt.mtPlus==None:
						el.nodes[1].elPlus = mc.elements[mt.mtPlus.e0]
						mc.nodes[mt.mtPlus.n0].elMinus = el
					else:
						el.nodes[1].elPlus = None;
					"""
					print "After polymerizing..."
					print mt
					for nid in range(mt.n0,mt.n2+3):
						print mc.nodes[nid]
					for eid in range(mt.e0,mt.e2+1):
						print mc.elements[eid]
					"""
		elif(mt.state==Depolymerizing): #Event may happen in this MT
			mc.mtRestore.append({'mt':mt,'n1':mt.n1,'e1':mt.e1, 'state':mt.state, \
								    'timeToNextEvent':mt.timeToNextEvent})

			for eid in range(mt.e1,mt.e0-1,-1): # Loop from + to - side 
				el = mc.elements[eid]
				if el.timeToNextEvent<=0:
					"""
					print "Before depolymerizing..."
					print mt
					for nid in range(mt.n0,mt.n2+3):
						print mc.nodes[nid]
					for eid in range(mt.e0,mt.e2+1):
						print mc.elements[eid]
					"""
					mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
					# Set new state and timeToNextEvent in this element
					el.state = MicrotubuleInactive
					el.timeToNextEvent = float("inf")

					# Update pointers
					if(mt.e1<mt.e2):
						el.nodes[1].elPlus = mc.elements[mt.e1+1]
					else:
						el.nodes[1].elPlus = None;

					if not mt.mtPlus==None:
						el.nodes[0].elPlus = mc.elements[mt.mtPlus.e0] 
						mc.nodes[mt.mtPlus.n0].elMinus = mc.elements[mt.e1-1]
					else:
						el.nodes[0].elPlus = None

					
					# Update n1 and e1 in this element
					mt.n1 -=1
					mt.e1 -=1	
					"""
					print "After depolymerizing..."
					print mt
					for nid in range(mt.n0,mt.n2+3):
						print mc.nodes[nid]
					for eid in range(mt.e0,mt.e2+1):
						print mc.elements[eid]
					"""

	cpdef ChangeMTState(self, MT mt, ModelContainer mc, Solver s):
		cdef Element el
		cdef int eid, nMTCS
		cdef double x1, xec, randMT,t_cycle,p_poly,p_stat

		mc.mtRestore.append({'mt':mt,'n1':mt.n1,'e1':mt.e1, 'state':mt.state, \
								    'timeToNextEvent':mt.timeToNextEvent})

		#print "Microtubule needs (de)polymerization."
		#print "Before: ", mt
		#for eid in range(mt.e0,mt.e2+1):
		#	print mc.elements[eid]

		"""
		# Compute new state and timeToNextEvent for this MT
		randMT = np.random.rand()
		t_cycle = self.tMTpoly+self.tMTstat+self.tMTdepoly
		p_poly = self.tMTpoly / t_cycle;
		p_stat = self.tMTstat /t_cycle;
		try:
			if np.abs(mc.currLMT-self.goalLMT)>1.e-4:
				p_poly = (1-self.fracLMT)*p_poly+self.fracLMT*(mc.currLMT<self.goalLMT)
				p_stat = (1-self.fracLMT)*p_stat;
		except:
			pass

		if randMT< p_poly: # Polymerizing
			mt.state = Polymerizing
			mt.timeToNextEvent = np.random.rand()*self.tMTpoly
		elif randMT<(p_poly+p_stat): # Stationary
			mt.state = Stationary
			mt.timeToNextEvent = np.random.rand()*self.tMTstat
		else: # Depolymerizing
			mt.state = Depolymerizing
			mt.timeToNextEvent = np.random.rand()*self.tMTdepoly
		"""
		
		"""
		# Compute new state and timeToNextEvent for this MT
		randMT = np.random.rand()
		if randMT<self.fracLMT and not mc.currLMT==self.goalLMT:		
			if mc.currLMT<self.goalLMT:
				mt.state = Polymerizing
				mt.timeToNextEvent = np.random.rand()*self.tMTpoly
			else:
				mt.state = Depolymerizing
				mt.timeToNextEvent = np.random.rand()*self.tMTdepoly
		"""
		if mt.state==Stationary: # Polymerizing
			# Compute number of MT in cross section at distal end of mt
			nMTCS = self.GetNumMTinCS(mt,mc,s)
			mt.state = Polymerizing
			mt.timeToNextEvent = (1.-self.fracLMT*(nMTCS/self.goalNMT-1.))\
			                     *np.random.rand()*self.tMTpoly
		elif mt.state == Depolymerizing: # Stationary
			mt.state = Stationary
			mt.timeToNextEvent = np.random.rand()*self.tMTstat
		else: # Depolymerizing
			mt.state = Depolymerizing
			mt.timeToNextEvent = np.random.rand()*self.tMTdepoly	
		

		# Change timeToNextEvent for all elements in this MT
		# Note, this is based on the reference configuration
		x1 = mc.nodes[mt.n1].x
		for eid in range(mt.e0,mt.e2+1):
			el = mc.elements[eid]
			mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
			
			xec = 0.5*(el.nodes[0].x+el.nodes[1].x) # Center of element
			
			if eid==mt.e0:	# First element never depolymerizes
				el.timeToNextEvent = float("inf")
			elif mt.state==Polymerizing:
				if el.state == MicrotubuleInactive:
					el.timeToNextEvent = (xec-x1)/self.MTpolyRate
				else:
					el.timeToNextEvent = float("inf")
			elif mt.state==Stationary:
				el.timeToNextEvent = float("inf")
			else:
				if el.state == MicrotubuleInactive or eid==mt.e0:
					el.timeToNextEvent = float("inf")
				else:
					el.timeToNextEvent = (x1-xec)/self.MTdepolyRate

		#print "After:  ", mt
		#for eid in range(mt.e0,mt.e2+1):
		#	print mc.elements[eid] 

	cpdef GetNumMTinCS(self, MT mt, ModelContainer mc, Solver s):
		cdef MT mti
		cdef list mtl
		cdef int nMTCS = 0
		cdef double x0,x1,xmt
		cdef Node n0,n1,nmt

		# Compute x location of node n1 of this mt
		nmt  = mc.nodes[mt.n1]
		xmt  = nmt.x+nmt.Dof(s.dc)[0]

		for mtl in mc.MTs:
			for mti in mtl:
				n0 = mc.nodes[mti.n0]
				n1 = mc.nodes[mti.n1]

				x0 = n0.x+n0.Dof(s.dc)[0]
				x1 = n1.x+n1.Dof(s.dc)[0]

				if xmt>x0 and xmt<x1:
					nMTCS+=1
					break

		return nMTCS