# -*- coding: utf-8 -*-
# In this mechanism, elements reattach to the same nodes they came from, i.e. no search.
# It is then check that the element does not connect to an inactive mt
from Mechanism cimport *
from CBARX cimport *
from ModelContainer cimport *
from Solver cimport *

import numpy as np
cimport numpy as np
from math import *

cdef class MECH_EL01(Mechanism):
	def __init__(self,tCrea, tDest):
		super(MECH_EL01,self).__init__();

		self.tCrea 			= tCrea
		self.tDest 			= tDest


	cpdef Apply(self, Element el, ModelContainer mc, Solver s):
		cdef int i
		cdef int nCreate, ndof, countK
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
			
			if  el.state==MotorAttachedSlack or el.state==MotorJustAttachedSlack:
				# Destroy element by changing its connectivity to the storage nodes of the model.
				el.dummyNodes = el.nodes
				el.nodes = mc.storageNodes
				el.timeToNextEvent = np.random.rand()*self.tCrea
				el.state           = MotorJustReleased
				nCreate -= 1
			elif el.state==MotorReleased or el.state==MotorJustReleased:
				# Create element by changing its connectivity to its dummyNodes
				el.nodes = el.dummyNodes
				el.dummyNodes = []
				el.timeToNextEvent = np.random.rand()*self.tDest
				el.state           = MotorJustAttachedSlack
				nCreate += 1
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