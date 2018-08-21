# -*- coding: utf-8 -*-
from Solver cimport *
from Element cimport *
from LOAD cimport *
from MPC cimport *
from SPC cimport *
cimport OutputHelper as oh
import numpy as np
cimport numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
from math import *
import copy, time
from multiprocessing import Array, Process

cdef class SOL2X(Solver):
	def __init__(self,mc,dc,doParallel,nProc=4):
		super(SOL2X,self).__init__('SOL2X',mc,dc,doParallel,nProc);
		self.tolerance = 1.e-6	# Default tolerance
		self.plotOutput = False

	cpdef Solve(self):
		cdef DataContainer dc = self.dc
		cdef ModelContainer mc = self.mc
		cdef Element el
		cdef Node nod
		cdef LOAD load
		cdef SPC s
		cdef int iteration, i, j, iel, iload, numdof, idK0, idK1, idR0, idR1, nElPerProc,nLoadPerProc
		cdef double residual0, residual, tolerance
		cdef bint reset, mpcPresent
		cdef list dof, dofs, dofm, dofr, dofp, dup
		cdef list rowK, colK, rowR, colR
		cdef list iel0, iel1, iload0, iload1, pList
		cdef np.ndarray[double, ndim=1] du

		# Initialize data
		for el in self.mc.elements:
			el.InitializeData(dc)
		"""
		for nod in self.mc.nodes:
			try:
				nod.InitializeCurvature(dc)
			except:
				pass
		"""
		# Share initialized data
		dc.ShareData()

		# Check whether we have MPC
		mpcPresent = mc.nmpc>0

		# Get id of dof sets (s=slave, m=master, r=free, p=prescribed)
		if mpcPresent:
			dofs, dofm, dofr, dofp, M = self.MPCMatrix(mc)
		else:
			dofp=[]
			for s in mc.spc:
				dofp.extend([s.node.dofID[s.dofIDs[i]] for i in range(len(s.dofIDs))])
			dofr = range(mc.ndof)
			dofr = np.delete(dofr,dofp).tolist()

		# Compute vector row, col for StiffnessMatrix and ResidualForceVector
		self.RowColumnVectors();
		
		# Shared data vectors
		#self.datK = np.zeros(len(self.rowK))
		#self.datR = np.zeros(len(self.rowR))
		self.datK = Array('d', len(self.rowK), lock=False)
		self.datR = Array('d', len(self.rowR), lock=False)

		if(self.doParallel):

			# Elements
			nElPerProc= mc.nel/self.nProc
			iel0=[0]
			iel1=[]
			for i in range(1,self.nProc):
				if i<=mc.nel%self.nProc:
					iel0.append(iel0[-1]+nElPerProc+1)
				else:
					iel0.append(iel0[-1]+nElPerProc)
				iel1.append(iel0[-1])
			iel1.append(mc.nel)

			# Loads
			nLoadPerProc= mc.nloads/self.nProc
			iload0=[0]
			iload1=[]
			for i in range(1,self.nProc):
				if i<=mc.nloads%self.nProc:
					iload0.append(iload0[-1]+nLoadPerProc+1)
				else:
					iload0.append(iload0[-1]+nLoadPerProc)
				iload1.append(iload0[-1])
			iload1.append(mc.nloads)
		
		# Initialize values
		dc.step = 0;
		dc.time = 0.;
		if self.plotOutput:
			oh.ParaviewOutput(mc,dc) # Write paraview output
		mc.WriteStepOutput(self) # Write model-specific output for this step

		dc.dt   = self.dt0;
		dc.dt0  = self.dt0;
		dc.time = dc.dt;

		while(dc.dt>1e-6 and dc.step<self.maxStep):
			dc.step+=1
			iteration = 0;
			residual = 1.;
			reset = False;

			# Decrement timeToNextEvent in each element and update model
			mc.UpdateModel(self)
			mc.RestoreModel(self)
			mc.UpdateModel(self)

			# Update dof0, stretch0 and Fg0, dt0
			for id1 in range(len(dc.Fg_)):
				dc.Fg0_[id1] = dc.Fg_[id1]

			for id1 in range(len(dc.stretch_)):
				dc.stretch0_[id1] = dc.stretch_[id1]

			for id1 in range(len(dc.curvature_)):
				dc.curvature0_[id1] = dc.curvature_[id1]

			for id1 in range(len(dc.theta_)):
				dc.theta0_[id1] = dc.theta_[id1]

			for id1 in range(mc.ndof):						
				dc.dof0[id1]=dc.dof[id1]

			for id1 in range(len(dc.Rp)):						
				dc.Rp0[id1]=dc.Rp[id1]
			
			dc.dt0 = dc.dt

			while(residual>self.tolerance and iteration< self.maxIter):
				iteration+= 1;

				# Stiffness matrix
				if(self.doParallel):
					pList = []
					[pList.append(Process(target=self.StiffnessData, args=(range(iel0[i],iel1[i]),))) for i in range(self.nProc)]
					[pList[i].start() for i in range(self.nProc)]
					[pList[i].join() for i in range(self.nProc)]
				else:
					self.StiffnessData(range(self.mc.nel))
				K = sps.coo_matrix((self.datK,(self.rowK,self.colK)),(mc.ndof,mc.ndof));
				if mpcPresent:
					Kmm = K.tocsr()[dofm,:].tocsc()[:,dofm]
					Kms = K.tocsr()[dofm,:].tocsc()[:,dofs]
					Kss = K.tocsr()[dofs,:].tocsc()[:,dofs]
					K = Kmm+Kms.dot(M)+M.T.dot(Kms.T)+M.T.dot(Kss.dot(M))

				# Global residual in first iteration
				if(iteration==1):
					if(self.doParallel):
						pList = []
						[pList.append(Process(target=self.ResidualDataExt, args=(range(iload0[i],iload1[i]),dc.time,))) for i in range(self.nProc)]
						[pList.append(Process(target=self.ResidualDataInt, args=(range(iel0[i],iel1[i]),))) for i in range(self.nProc)]
						[pList[i].start() for i in range(2*self.nProc)]
						[pList[i].join() for i in range(2*self.nProc)]
					else:
						# External loads
						self.ResidualDataExt(range(mc.nloads),dc.time)
						# Minus internal loads
						self.ResidualDataInt(range(mc.nel))

					R = sps.coo_matrix((self.datR,(self.rowR,self.colR)),(mc.ndof,1));
					if mpcPresent:
						Rm  = R.tocsr()[dofm,:]
						Rs  = R.tocsr()[dofs,:]
						R = Rm + M.T.dot(Rs)

				# Apply SPC
				dup  = [];
				for s in mc.spc:
					value = s.Get(dc.time,dc.dt)

					if(iteration==1):
						dup.extend([value for i in range(len(s.dofIDs))])
					else:
						dup.extend(np.zeros(len(s.dofIDs)))

				# Slice and solve
				Krr = K.tocsr()[dofr,:].tocsc()[:,dofr]
				Krp = K.tocsr()[dofr,:].tocsc()[:,dofp]

				Rr  = R.tocsr()[dofr,:]
				RHS = Rr-Krp.dot(np.matrix(dup).T)
				dur = spsl.spsolve(Krr,RHS,);

				du = np.zeros(mc.ndof);
				if mpcPresent:
					dum = np.zeros(mc.ndof-mc.nmpc)
					dum[dofr] = dur
					dum[dofp] = dup
					du[dofm] = dum
					du[dofs] = M.dot(np.matrix(dum).T)#.flatten()
				else:
					du[dofr] = dur
					du[dofp] = dup
				dc.dof = dc.dof + du

				"""
				print "Total left wall force  = ", np.sum(np.array(dc.Rp)[0:10])
				print "Total stor.nodes force = ", np.sum(np.array(dc.Rp)[10:12])
				print "Total visc.nodes force = ", np.sum(np.array(dc.Rp)[12:31])
				print "Total reaction force   = ", np.sum(np.array(dc.Rp)[0:10])+np.sum(np.array(dc.Rp)[12:31])
				print "Total pulling force    = ", np.array(dc.Rp)[31]
				print "Total force error      = ", np.sum(np.array(dc.Rp)[0:10])+np.sum(np.array(dc.Rp)[12:32])
				"""
				
				# Assemble global residual
				if(self.doParallel):
					pList = []
					[pList.append(Process(target=self.ResidualDataExt, args=(range(iload0[i],iload1[i]),dc.time,))) for i in range(self.nProc)]
					[pList.append(Process(target=self.ResidualDataInt, args=(range(iel0[i],iel1[i]),))) for i in range(self.nProc)]
					[pList[i].start() for i in range(2*self.nProc)]
					[pList[i].join() for i in range(2*self.nProc)]
				else:
					# External loads
					self.ResidualDataExt(range(mc.nloads),dc.time)
					# Minus internal loads
					self.ResidualDataInt(range(mc.nel))
				R = sps.coo_matrix((self.datR,(self.rowR,self.colR)),(mc.ndof,1));
				if mpcPresent:
						Rm  = R.tocsr()[dofm,:]
						Rs  = R.tocsr()[dofs,:]
						R = Rm + M.T.dot(Rs)
				dc.Rp  = R.tocsr()[dofp,:].toarray().flatten()
				
				residual0 = residual
				residual = sqrt((R.tocsr()[dofr,:].data ** 2).sum())/mc.ndof

				oh.WriteToLog(mc,"Step: "+str(dc.step)+"  \t Iter: "+str(iteration)+"   \t dt: "+'{:6.4f}'.format(dc.dt)+"   \t t: "+'{:6.4f}'.format(dc.time)+"   \t Residual: "+'{:6.4e}'.format(residual)+"   \t SimDef: "+'{:6.4e}'.format(max(np.array([dc.dof[i] for i in range(1,mc.ndof,2)]))))
				#if (residual > residual0):
				#	oh.WriteToLog(mc,"Residual went up, time step decreased.\n")
				#	reset = True;
				#	break;	

				"""	
				np.set_printoptions(suppress=True)
				#print "\n\nKrr          = \n", Krr.todense()
				print "Rr           = \n", R.tocsr()[dofr,:]
				print "dur          = \n", dur
				"""

				if np.isnan(residual):
					raise("ERROR: Residual isnan!!")

			# Compute new time step
			if(iteration==self.maxIter and iteration!=1 and residual > self.tolerance):
				oh.WriteToLog(mc,"Too many iterations, time step decreased.")
				reset = True

			if(reset):
				dc.step-=1
				dc.time-=dc.dt 

				# Restore dof, Rp, stretch and Fg, dt
				for id1 in range(mc.ndof):						
					dc.dof[id1]=dc.dof0[id1]

				for id1 in range(len(dc.Rp)):						
					dc.Rp[id1]=dc.Rp0[id1]

				for id1 in range(len(dc.Fg_)):
					dc.Fg_[id1] = dc.Fg0_[id1]

				for id1 in range(len(dc.stretch_)):
					dc.stretch_[id1] = dc.stretch0_[id1]

				for id1 in range(len(dc.curvature_)):
					dc.curvature_[id1] = dc.curvature0_[id1]

				for id1 in range(len(dc.theta_)):
					dc.theta_[id1] = dc.theta0_[id1]

				dc.dt = dc.dt0

				# Restore the model and timeToNextEvent in each element
				mc.RestoreModel(self)

				dc.dt/=2.
				dc.time+=dc.dt 

				# Check if time increment becomes too small
				if(dc.dt<=self.dtMin):
					oh.WriteToLog(mc,"ERROR: Time increment too small.\n")
					return
			else:
				if self.plotOutput:
					oh.ParaviewOutput(mc,dc) # Write paraview output
				mc.WriteStepOutput(self) # Write model-specific output for this step

				if(iteration<=self.maxIterInc):
					dc.dt*=1.25
				if(dc.dt>=self.dtMax):
					dc.dt = self.dtMax
				if(dc.time+dc.dt>self.tEnd):
					dc.dt = self.tEnd-dc.time

				dc.time+=dc.dt 

				oh.WriteToLog(mc,"\n")


	cpdef StiffnessData(self, list eids):
		cdef int i,j
		cdef np.ndarray[double, ndim=1] kel
		cdef Element el
		
		for j in eids:
			el = self.mc.elements[j]
			kel = el.BuildElementMatrix(self.dc).flatten('F')
			for i in range(el.datKID[0],el.datKID[1]):
				self.datK[i] = kel[i-el.datKID[0]]

	cpdef ResidualDataExt(self, list lids, double time):
		cdef int i,j
		cdef np.ndarray[double, ndim=1] rel
		cdef LOAD load

		# External loads
		for j in lids:
			load = self.mc.loads[j]
			rel = load.Get(time)
			for i in range(load.datRID[0],load.datRID[1]):
				self.datR[i] = rel[i-load.datRID[0]]

	cpdef ResidualDataInt(self, list eids):
		cdef int i,j
		cdef np.ndarray[double, ndim=1] rel
		cdef Element el
		
		# Minus internal loads
		for j in eids:
			el = self.mc.elements[j]
			rel = -1.0*el.BuildInternalForceVector(self.dc)
			for i in range(el.datRID[0],el.datRID[1]):
				self.datR[i] = rel[i-el.datRID[0]]

	cpdef MPCMatrix(self, ModelContainer mc):
		cdef list dofs, dofm, dofr, dofp
		cdef list rowM, colM, datM
		cdef MPC m
		cdef SPC s

		dofs = []; dofm = []; dofr = []; dofp = [];
		rowM = []; colM = []; datM = [];

		# get slave dof, master dof, and weight matrix M
		count = 0;
		for m in mc.mpc:
			dofs.append(m.slaveDofID)
			rowM.extend([count for i in range(len(m.masterDofIDs))])
			colM.extend(m.masterDofIDs);
			datM.extend(m.masterWeights);
			count = count+1;

		# dofm is all dof except for dofs
		dofm = range(mc.ndof)
		dofm = np.delete(dofm,dofs).tolist()

		# M such that u_s = M*u:
		M = sps.coo_matrix((datM,(rowM,colM)),(count,mc.ndof));
		# M such that u_s = M*u_m:
		M = M.tocsc()[:,dofm]

		# Compute free and prescribed dof indices dofr and dofp within dofm
		for s in mc.spc:
			try:
				dofp.extend([dofm.index(s.node.dofID[s.dofIDs[i]]) for i in range(len(s.dofIDs))])
			except:
				raise("ERROR: SPC cannot be applied to slave node of MPC!!")
		dofr = range(mc.ndof-mc.nmpc)
		dofr = np.delete(dofr,dofp).tolist()

		return dofs, dofm, dofr, dofp, M

	cpdef RowColumnVectors(self):
		cdef Element el 
		cdef LOAD load
		cdef int idK0 = 0, idK1 = 0, idR0=0, idR1=0
		cdef int i, j

		self.rowK=[]; self.colK=[]; self.rowR = [];

		for el in self.mc.elements:
			dof = el.DofID().tolist();

			# Stiffness matrix
			self.rowK.extend([i for j in range(len(dof)) for i in dof])
			self.colK.extend([i for i in dof for j in range(len(dof))])
			idK1 = idK0+len(dof)*len(dof);
			el.datKID = [idK0, idK1]
			idK0 = idK1;

			# Residual vector - internal load
			self.rowR.extend(dof)
			idR1 = idR0+len(dof);
			el.datRID = [idR0, idR1]
			idR0 = idR1;

		for load in self.mc.loads:
			# Residual vector - external load
			dof = load.node.dofID;
			self.rowR.extend(dof)
			idR1 = idR0+len(dof);
			load.datRID = [idR0, idR1]
			idR0 = idR1;
		self.colR = np.zeros(len(self.rowR)).tolist()
