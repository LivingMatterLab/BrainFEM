from ModelContainer cimport *
from Amplitude cimport *
from MAT_NEOH cimport * 
from MAT_LE cimport *
from PBAR cimport *
from PBAR1 cimport *
from Node cimport *
from CBAR cimport *
from CBAR1 cimport *
from SPC cimport *
from LOAD cimport *

cimport OutputHelper as oh
cimport ElementHelper as eh

from math import *
import numpy as np
cimport numpy as np

cdef class BarTest_B3(ModelContainer):
	def __init__(self):
		super().__init__()
		print "Initializing BarTest_B3"

	cpdef BuildModel(self, object p):
		self.numDofPerNode = 3

		# INPUT
		# X, Y, Z location of first node
		cdef double Xn = 1000.0;
		cdef double Yn = 0.0;
		cdef double Zn = 0.0;
		cdef double Lbar = sqrt(Xn**2+Yn**2+Zn**2)
		cdef double Lrest = 1.*Lbar;					# rest length

		# Deformation of second node
		#cdef double Vn = 100.;
		#cdef double Un = 100.;
		#cdef double Wn = -1000.;

		# Force/displacement at second node
		cdef double Fx = 75.;
		cdef double Vn = Xn/1.;
		cdef double Wn = Xn/1.;
		cdef double Fmotor = 0.;			# Motor force in element

		# Create amplitude
		cdef double tEnd = 15.;
		cdef double tLoad = 10.;
		self.amplitudes.append(Amplitude(np.array([0.,tLoad,tEnd]),np.array([0.,1.,1.])));

		# Create a material
		cdef dict matProp = {'E':10.0,'nu':0.4}
		self.materials.append(MAT_LE(matProp))
		self.materials[0].localID  = 12;
		print self.materials[0]

		# Create a PBAR property
		cdef dict propProp = {'area':5.,'force':Fmotor}
		self.properties.append(PBAR1(self.materials[0],propProp,self.amplitudes[0]))
		self.properties[0].localID = 4;
		print self.properties[0]

		# Create 2 nodes
		self.nodes.append(Node([0.,0.,0.]));
		self.nodes[0].localID = 0;
		self.nodes[0].dofID = range(self.numDofPerNode);
		self.nodes.append(Node([Xn,Yn,Zn]));
		self.nodes[1].localID = 1;
		self.nodes[1].dofID = range(self.numDofPerNode,2*self.numDofPerNode);

		# Create CBAR element
		self.elements.append(CBAR([self.nodes[0], self.nodes[1]],self.properties[0],Lrest));
		self.elements[0].localID = 0;
		print self.elements[0]

		# SPC
		# Clamp node 0
		self.spc.append(SPC(self.nodes[0],range(self.numDofPerNode),0.,self.amplitudes[0]));
		# Move node 1
		#self.spc.append(SPC(self.nodes[1],[0],Un));
		#self.spc.append(SPC(self.nodes[1],[1],Vn));
		#self.spc.append(SPC(self.nodes[1],[2],Wn));

		# Apply load to node 1, and constrain y,z directions
		self.loads.append(LOAD(self.nodes[1],Fx*np.array([1.,0.,0.]),self.amplitudes[0]));
		self.spc.append(SPC(self.nodes[1],[1],Vn,self.amplitudes[0]));
		self.spc.append(SPC(self.nodes[1],[2],Wn,self.amplitudes[0]));

		# length of all matrices
		self.SetListLengths()

	cpdef TestElementMatrix(self,DataContainer dc):
		# Here we test whether the stiffness matrix has stiffness orthogonal to 
		# bar element in presence of internal stresses
		print "Test Element Matrix (works when node[1] is on X-axis)"

		# Set stretch
		cdef double stretchTest = 1.5;
		cdef double L = np.linalg.norm(np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc))
		cdef double area = self.properties[0].area
		cdef double mu  = self.materials[0].mu
		dc.stretch_ = stretchTest*np.ones(1)
		dc.stretch0_ = stretchTest*np.ones(1)
		dc.dt = 1.
		self.elements[0].stretchID = [0]
		self.elements[0].stretch0ID = [0]
		dc.dof[3] = (stretchTest-1.)*L

		# Current method
		kel = self.elements[0].BuildElementMatrix(dc)
		print "CBAR: kel = \n" , kel

		# Use second piola stress, and method from delft thesis
		dNmat = np.array([[1.,0.,0.,-1.,0.,0.],[0.,1.,0.,0.,-1.,0.],[0.,0.,1.,0.,0.,-1.]])/L
		Pmat = dNmat.T.dot(dNmat)


		rho0 = np.array([0.,0.,0.,L,0.,0.])
		v    = np.array([0.,0.,0.,L*(stretchTest-1.),0.,0.])

		Bmat = (rho0+v).T.dot(Pmat)

		Me = self.materials[0].Piola1Stiffness_1d(stretchTest)
		s11 = mu*(1.-1./stretchTest**3)   # 2nd piola stress
		c11 = 3*mu/stretchTest**5	      # ds11/de11, e11 = 0.5*(stretchTest**2-1)

		mat2 = area*L*(s11*Pmat+c11*np.outer(Bmat,Bmat))

		print "Delft Thesis method: kel = \n" , mat2

		# Compute error
		errMat = np.linalg.norm(kel-mat2)
		print "Error =  ", errMat


	cpdef PostProcess(self, DataContainer dc):
		# Repeat input:
		# X, Y, Z location of first node
		cdef double Xn = 1000.0;
		cdef double Yn = 0.0;
		cdef double Zn = 0.0;
		cdef double Lbar = sqrt(Xn**2+Yn**2+Zn**2)

		# Force at second node
		cdef double Fx = 75;
		cdef double Vn = Xn/1.;
		cdef double Wn = Xn/1.;

		# Compute stretch
		ex0 = np.array(self.nodes[1].loc)-np.array(self.nodes[0].loc);
		ex  = ex0 + np.array(self.nodes[1].Dof(dc))-np.array(self.nodes[0].Dof(dc))
		l   = np.linalg.norm(ex);
		stretch = l/self.elements[0].restLength;
		stretch_rate = 0.0;

		# Compute force
		M = self.properties[0].Piola1Stiffness(stretch,stretch_rate,dc);
		F = M['P']
		angle = np.arctan(np.sqrt(Vn**2+Wn**2)/ex[0])
		Fxc = F*np.cos(angle)
		
		oh.WriteToLog(self,'\n PostProcess:')
		oh.WriteToLog(self,'ex = ' + str(ex))
		oh.WriteToLog(self,'stretch = ' + str(stretch))
		oh.WriteToLog(self,'F = ' + str(F))
		oh.WriteToLog(self,'Fxc = ' + str(Fxc))
		oh.WriteToLog(self,'error = ' + str((Fxc-Fx)/Fx))


