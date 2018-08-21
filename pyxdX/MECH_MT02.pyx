# Mechanism that allows MT to polymerize, depolymerize and to be stationary
# Additionally, this MT will be 'attached' by stiff crosslink to growth cone
# if it becomes to close to it.
from MECH_MT01 cimport *
from ModelContainer cimport *
from Element cimport *
from CBARX2 cimport *
from Solver cimport *

import numpy as np
cimport numpy as np
from math import *
cimport ElementHelper as eh

cdef class MECH_MT02(MECH_MT01):
	def __init__(self,tMTpoly,tMTstat,tMTdepoly,MTpolyRate,MTdepolyRate,thresDist,stabDistalMT,fracLMT=0.,goalNMT=9.5):
		super(MECH_MT02,self).__init__(tMTpoly,tMTstat,tMTdepoly,MTpolyRate,MTdepolyRate,fracLMT,goalNMT);

		self.stabDistalMT = stabDistalMT
		self.thresDist    = thresDist

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
		if mt.mtMinus==None and mc.nodes[mt.n0].x<1.e-6:
			# First MT, if attached to wall, does never polymerize or depolymerize
			mt.state = Stationary
			mt.timeToNextEvent = float("inf")
		if self.stabDistalMT and mt.mtPlus==None and mc.nodes[mt.n1].x-1.e-6>mc.lAxon:
			# Last MT, if reaching to end of axon and stabDistalMT, does never polymerize or depolymerize
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
		
		# Initialize GC crosslink element
		self.InitializeGCNode(mt, mc)

	cpdef Apply(self, MT mt, ModelContainer mc, Solver s):
		cdef Element el
		cdef int eid
		cdef double x1, xec, randMT

		# Apply GC crosslink
		self.ApplyGCNode(mt,mc,s)

		if(mt.timeToNextEvent<=0):
			self.ChangeMTState(mt,mc,s)
			self.ApplyGCNode(mt,mc,s)
		elif(mt.state==Stationary): # No event is happening in this MT
			pass;
		elif(mt.state==Polymerizing): #Event may happen in this MT
			mc.mtRestore.append({'mt':mt,'n1':mt.n1,'e1':mt.e1, 'state':mt.state, \
								    'timeToNextEvent':mt.timeToNextEvent})
			
			for eid in range(mt.e1,mt.e2+1): # Loop from - to + side 
				el = mc.elements[eid]
				if el.timeToNextEvent<=0:
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

					# Check GC Node
					self.ApplyGCNode(mt,mc,s)

		elif(mt.state==Depolymerizing): #Event may happen in this MT
			mc.mtRestore.append({'mt':mt,'n1':mt.n1,'e1':mt.e1, 'state':mt.state, \
								    'timeToNextEvent':mt.timeToNextEvent})

			for eid in range(mt.e1,mt.e0-1,-1): # Loop from + to - side 
				el = mc.elements[eid]
				if el.timeToNextEvent<=0:
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

		# Apply GC crosslink element
		self.ApplyGCNode(mt,mc,s)
		
		return 0

	cpdef InitializeGCNode(self, MT mt, ModelContainer mc):
		cdef Element el
		cdef double distNodes

		# Initialize GC crosslink element	
		if mt.eGC > -1:
			el = mc.elements[mt.eGC]
			el.restLength = np.sqrt( (el.nodes[0].y-el.nodes[1].y)**2 + \
								     (el.nodes[0].z-el.nodes[1].z)**2 + \
								     self.thresDist**2)

			# Check if elment needs to be attached or detached
			distNodes = np.abs( el.nodes[0].x-el.nodes[1].x)
			if(distNodes<self.thresDist-1.e-6):
				# Attached
				el.state = MotorAttachedSlack	
			else:
				el.nodes = mc.storageNodes
				el.state = MotorReleased	
			el.timeToNextEvent = 1.e9
			
	cpdef ApplyGCNode(self, MT mt, ModelContainer mc, Solver s):
		cdef Element el,eli
		cdef double distNodes, elStretch
		cdef list nodDof0, nodDof1

		# Check GC crosslink element
		if mt.eGC > -1:
			el = mc.elements[mt.eGC]
			if el.state==MotorJustAttachedSlack:
				mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
				el.state = MotorAttachedSlack;
			if el.state==MotorJustReleased:
				mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent})
				el.state = MotorReleased;

			# Compute distance of element
			nodDof0 = mc.nodes[mt.n1].Dof(s.dc)
			nodDof1 = mc.growthConeNode.Dof(s.dc)
			distNodes = np.abs( mc.growthConeNode.x+nodDof1[0] \
							   -mc.nodes[mt.n1].x-nodDof0[0])

			# Compute element stretch
			elStretch = el.CurrentLength(s.dc)/el.restLength

			#print distNodes, self.thresDist
			# Attach crosslink if distance below threshold
			if(distNodes<self.thresDist-1.e-6 and \
				(el.state==MotorJustReleased or el.state==MotorReleased)):
				mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
				                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
				
				# Attach element by changing its connectivity to its dummyNodes
				el.nodes = [mc.nodes[mt.n1], mc.growthConeNode]
				el.timeToNextEvent = 1.e9
				el.state           = MotorJustAttachedSlack

				# Update rowR, rowK, colK in the solver
				eh.UpdateConnectivities(el,s)
			
			elif((distNodes>1.001*self.thresDist or elStretch>1.) and \
				(el.state==MotorJustAttachedSlack or el.state==MotorAttachedSlack)):
				mc.elRestore.append({'el':el,'state':el.state,'timeToNextEvent':el.timeToNextEvent,'restLength':el.restLength,\
				                     'nodes':[n.localID for n in el.nodes] ,'dummyNodes':[n.localID for n in el.dummyNodes]})
				
				# Detach element by changing its connectivity to its dummyNodes
				el.nodes = mc.storageNodes
				el.timeToNextEvent = 1.e9
				el.state           = MotorJustReleased

				# Update rowR, rowK, colK in the solver
				eh.UpdateConnectivities(el,s)
			
			# If Microtuble was polymerizing, then prevent 
			# activation of new elements that
			if (mt.state == Polymerizing and \
				(el.state==MotorJustAttachedSlack or el.state==MotorAttachedSlack)):
				for eid in range(mt.e1,mt.e2+1): # Loop from - to + side 
					eli = mc.elements[eid]
					mc.elRestore.append({'el':eli,'state':eli.state,'timeToNextEvent':eli.timeToNextEvent})
					eli.timeToNextEvent = float("inf")
