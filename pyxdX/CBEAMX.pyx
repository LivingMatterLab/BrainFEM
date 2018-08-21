# -*- coding: utf-8 -*-
from CBEAM cimport *
cimport numpy as np
import numpy as np

cdef class CBEAMX(CBEAM):
	def __init__(self,nodes,elementProperty,state=-1,timeToNextEvent=1.e9, restLength=0.0):
		super(CBEAMX,self).__init__(nodes,elementProperty,restLength);
		self.state = state
		self.timeToNextEvent = timeToNextEvent;
		self.dummyNodes = []
		self.mechanism = None
		self.minusID = -1;			# Node on minus side (indicates direction of polarity)
		self.plusID = -1;			# Node on plus side (indicates direction of polarity)


	def __str__(self):
		cdef str nidD
		nidD = "None" if self.dummyNodes==[] else "[" + ', '.join(str(nod.localID) for nod in self.dummyNodes) + "]" 
		return super(CBEAMX,self).__str__() + "\t State = " + str(self.state) + "\t PlusID = " + str(self.plusID) + "\t TimeToNextEvent = " + str(self.timeToNextEvent)+ "\t dummyNodes = " + nidD

	cpdef double CurrentLength(self, DataContainer dc):
		cdef np.ndarray[double, ndim=1] ex
		cdef list nodDof0, nodDof1

		# Get nodal dof. Append 0's if necessary, i.e. if len(nod.dof)<len(nod.loc)
		nodDof0 = self.nodes[0].Dof(dc)
		while (len(nodDof0)<len(self.nodes[0].loc)): nodDof0.append(0.0) 

		nodDof1 = self.nodes[1].Dof(dc)
		while (len(nodDof1)<len(self.nodes[1].loc)): nodDof1.append(0.0) 

		ex = np.array(self.nodes[1].loc)     - np.array(self.nodes[0].loc) + \
		     np.array(nodDof1) - np.array(nodDof0)
		return np.linalg.norm(ex)