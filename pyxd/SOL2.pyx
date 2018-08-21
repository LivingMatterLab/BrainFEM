# -*- coding: utf-8 -*-
from Solver cimport *
from Element cimport *
from LOAD cimport *
from MPC cimport *
from SPC cimport *
from CONTACT4 cimport *
from PRESSURE cimport *
cimport OutputHelper as oh
import numpy as np
cimport numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
from math import *
import copy, time
from multiprocessing import Array, Process

cdef class SOL2(Solver):
	def __init__(self,mc,dc,doParallel,nProc=4):
		super(SOL2,self).__init__('SOL2',mc,dc,doParallel,nProc);
		self.plotOutput = False

	cpdef Solve(self):
		cdef DataContainer dc = self.dc
		cdef ModelContainer mc = self.mc
		cdef Element el
		cdef LOAD load
		cdef SPC s
		cdef int iteration, i, j, iel, iload, numdof, idK0, idK1, idR0, idR1, nElPerProc,nLoadPerProc,nContPerProc
		cdef double residual0, residual
		cdef bint reset, mpcPresent
		cdef list dof, dofs, dofm, dofr, dofp, dup
		cdef list rowK, colK, rowR, colR
		cdef list iel0, iel1, iload0, iload1, icont0, icont1, pList
		cdef np.ndarray[double, ndim=1] du


		# Initialize data
		for el in self.mc.elements:
			el.InitializeData(dc)

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
		rowK, colK, rowR, colR = self.RowColumnVectors();
		
		# Shared data vectors
		#self.datK = np.zeros(len(rowK))
		#self.datR = np.zeros(len(rowR))
		self.datK = Array('d', len(rowK), lock=False)
		self.datR = Array('d', len(rowR), lock=False)

		# Divide elements and loads over processes
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
			print iel0
			print iel1

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
			print iload0
			print iload1

			# Contacts
			nContPerProc= mc.ncontacts/self.nProc
			icont0=[0]
			icont1=[]
			for i in range(1,self.nProc):
				if i<=mc.ncontacts%self.nProc:
					icont0.append(icont0[-1]+nContPerProc+1)
				else:
					icont0.append(icont0[-1]+nContPerProc)
				icont1.append(icont0[-1])
			icont1.append(mc.ncontacts)
			print icont0
			print icont1

			# Pressures
			nPressPerProc= mc.npressures/self.nProc
			ipress0=[0]
			ipress1=[]
			for i in range(1,self.nProc):
				if i<=mc.npressures%self.nProc:
					ipress0.append(ipress0[-1]+nPressPerProc+1)
				else:
					ipress0.append(ipress0[-1]+nPressPerProc)
				ipress1.append(ipress0[-1])
			ipress1.append(mc.npressures)
			print ipress0
			print ipress1
		
		# Initialize values
		dc.step = 0;
		dc.time = 0.;
		if self.plotOutput:
			#oh.PlotOutput(mc, dc)	# Plot undeformed configuration
			oh.ParaviewOutput(self.mc,self.dc,0) # Write paraview output
		mc.WriteStepOutput(self) # Write model-specific output for this step
		dc.dt   = self.dt0;
		dc.time = dc.dt;
		
		while(dc.dt>1e-6 and dc.step<self.maxStep):
			dc.step+=1
			iteration = 0;
			residual = 1.;
			reset = False;

			# Update dof0, stretch0 and Fg0
			for id1 in range(len(dc.Fg_)):
				dc.Fg0_[id1] = dc.Fg_[id1]

			for id1 in range(len(dc.stretch_)):
				dc.stretch0_[id1] = dc.stretch_[id1]

			for id1 in range(len(dc.curvature_)):
				dc.curvature0_[id1] = dc.curvature_[id1]

			for id1 in range(len(dc.theta_)):
				dc.theta0_[id1] = dc.theta_[id1]

			for id1 in range(len(dc.Rp)):						
				dc.Rp0[id1]=dc.Rp[id1]

			for id1 in range(mc.ndof):						
				dc.dof0[id1]=dc.dof[id1]

			while(residual>1e-6 and iteration< self.maxIter):
				iteration+= 1;

				# Global residual in first iteration
				if(iteration==1):
					if(self.doParallel):
						pList = []
						[pList.append(Process(target=self.ResidualDataExt, args=(range(iload0[i],iload1[i]),range(ipress0[i],ipress1[i]),dc.time,))) for i in range(self.nProc)]
						#[pList.append(Process(target=self.ResidualDataCont, args=(range(icont0[i],icont1[i]),iteration,))) for i in range(self.nProc)]
						[pList.append(Process(target=self.ResidualDataInt, args=(range(iel0[i],iel1[i]),))) for i in range(self.nProc)]
						[pList[i].start() for i in range(2*self.nProc)]
						[pList[i].join() for i in range(2*self.nProc)]

						# Contact never parallel due to slaveIsContact vector not saved properly
						self.ResidualDataCont(range(mc.ncontacts),iteration)
						
					else:
						# External loads
						self.ResidualDataExt(range(mc.nloads),range(mc.npressures),dc.time)
						# Contact loads
						self.ResidualDataCont(range(mc.ncontacts),iteration)
						# Minus internal loads
						self.ResidualDataInt(range(mc.nel))
					R = sps.coo_matrix((self.datR,(rowR,colR)),(mc.ndof,1));
					if mpcPresent:
						Rm  = R.tocsr()[dofm,:]
						Rs  = R.tocsr()[dofs,:]
						R = Rm + M.T.dot(Rs)

				# Stiffness matrix
				if(self.doParallel):
					pList = []
					[pList.append(Process(target=self.StiffnessData, args=(range(iel0[i],iel1[i]),range(ipress0[i],ipress1[i]),))) for i in range(self.nProc)]
					#[pList.append(Process(target=self.StiffnessDataCont, args=(range(icont0[i],icont1[i]),))) for i in range(self.nProc)]
					[pList[i].start() for i in range(self.nProc)]
					[pList[i].join() for i in range(self.nProc)]

					# Contact never parallel due to slaveIsContact vector not shared properly
					self.StiffnessDataCont(range(self.mc.ncontacts))
				else:
					self.StiffnessData(range(self.mc.nel),range(self.mc.npressures))
					self.StiffnessDataCont(range(self.mc.ncontacts))
				K = sps.coo_matrix((self.datK,(rowK,colK)),(mc.ndof,mc.ndof));
				if mpcPresent:
					Kmm = K.tocsr()[dofm,:].tocsc()[:,dofm]
					Kms = K.tocsr()[dofm,:].tocsc()[:,dofs]
					Kss = K.tocsr()[dofs,:].tocsc()[:,dofs]
					K = Kmm+Kms.dot(M)+M.T.dot(Kms.T)+M.T.dot(Kss.dot(M))

				

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
				print 'K: '
				idof = [1, 7, 13,19,25,31,37,43,49,55,61,67]
				print K.tocsr()[idof,:].tocsc()[:,idof].todense()
				
				print 'Krr:'
				print Krr.todense()
				
				print '\n Rr:'
				print Rr.todense()
				print '\n Krp:'
				print Krp.todense()
				print '\n dur:'
				print dur
				print '\n dc.dof:'
				print np.asarray(dc.dof)
				"""

				
				
				# Assemble global residual
				if(self.doParallel):
					pList = []
					[pList.append(Process(target=self.ResidualDataExt, args=(range(iload0[i],iload1[i]),range(ipress0[i],ipress1[i]),dc.time,))) for i in range(self.nProc)]
					#[pList.append(Process(target=self.ResidualDataCont, args=(range(icont0[i],icont1[i]),iteration,))) for i in range(self.nProc)]
					[pList.append(Process(target=self.ResidualDataInt, args=(range(iel0[i],iel1[i]),))) for i in range(self.nProc)]
					[pList[i].start() for i in range(2*self.nProc)]
					[pList[i].join() for i in range(2*self.nProc)]

					# Contact never parallel due to slaveIsContact vector not shared properly
					self.ResidualDataCont(range(mc.ncontacts),iteration)
				else:
					# External loads
					self.ResidualDataExt(range(mc.nloads),range(mc.npressures),dc.time)
					# Contact loads
					self.ResidualDataCont(range(mc.ncontacts),iteration)
					# Minus internal loads
					self.ResidualDataInt(range(mc.nel))
				R = sps.coo_matrix((self.datR,(rowR,colR)),(mc.ndof,1));
				if mpcPresent:
						Rm  = R.tocsr()[dofm,:]
						Rs  = R.tocsr()[dofs,:]
						R = Rm + M.T.dot(Rs)
				dc.Rp  = R.tocsr()[dofp,:].toarray().flatten()
				
				residual0 = residual
				residual = sqrt((R.tocsr()[dofr,:].data ** 2).sum())/mc.ndof
				oh.WriteToLog(mc,"Step: "+str(dc.step)+"  \t Iter: "+str(iteration)+"   \t dt: "+'{:6.4f}'.format(dc.dt)+"   \t t: "+'{:6.4f}'.format(dc.time)+"   \t Residual: "+'{:6.4e}'.format(residual)+"   \t SimDef: "+'{:6.4e}'.format(max(np.array([dc.dof[i] for i in range(1,mc.ndof,3)]))))
				#if (residual > residual0):
				#	oh.WriteToLog(mc,"Residual went up, time step decreased.\n")
				#	reset = True;
				#	break;	

				if np.isnan(residual) and iteration==1:
					raise("ERROR: Residual isnan!!")
				elif np.isnan(residual):
					iteration = self.maxIter;
					oh.WriteToLog(mc,"Residual is NAN.\n")
					reset = True


				#if self.plotOutput:
				#	#oh.PlotOutput(mc, dc)	# Plot undeformed configuration
				#	oh.ParaviewOutput(self.mc,self.dc,iteration) # Write paraview output		

			# Compute new time step
			if(iteration==self.maxIter and iteration!=1 and residual > 1e-6):
				oh.WriteToLog(mc,"Too many iterations, time step decreased.\n")
				reset = True

			if(reset):
				dc.step-=1
				dc.time-=dc.dt 
				dc.dt/=2.
				dc.time+=dc.dt 
				# Restore dof0, stretch0 and Fg0
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

				# Check if time increment becomes too small
				if(dc.dt<=self.dtMin):
					oh.WriteToLog(mc,"ERROR: Time increment too small.\n")
					return
			else:
				if self.plotOutput:
					#oh.PlotOutput(self.mc,self.dc)	# Plot deformed configuration
					oh.ParaviewOutput(self.mc,self.dc,iteration) # Write paraview output
				mc.WriteStepOutput(self) # Write model-specific output for this step

				if(iteration<=self.maxIterInc):
					dc.dt*=1.25
				if(dc.dt>=self.dtMax):
					dc.dt = self.dtMax
				if(dc.time+dc.dt>self.tEnd):
					dc.dt = self.tEnd-dc.time

				dc.time+=dc.dt 

				oh.WriteToLog(mc,"\n")


	cpdef StiffnessData(self, list eids, list pids):
		cdef int i,j
		cdef np.ndarray[double, ndim=1] kel
		cdef Element el
		
		for j in eids:
			el = self.mc.elements[j]
			kel = el.BuildElementMatrix(self.dc).flatten('F')
			for i in range(el.datKID[0],el.datKID[1]):
				self.datK[i] = kel[i-el.datKID[0]]

		# Pressures
		for j in pids:
			pressure = self.mc.pressures[j]
			kel = pressure.BuildMatrix(self.dc).flatten('F')
			for i in range(pressure.datKID[0],pressure.datKID[1]):
				self.datK[i] = kel[i-pressure.datKID[0]]


	cpdef StiffnessDataCont(self, list cids):
		cdef int i,j
		cdef np.ndarray[double, ndim=1] kcont
		cdef CONTACT contact
		
		for j in cids:
			contact = self.mc.contacts[j]
			kcont = contact.BuildMatrix(self.dc).flatten('F')
			for i in range(contact.datKID[0],contact.datKID[1]):
				self.datK[i] = kcont[i-contact.datKID[0]]

	cpdef ResidualDataExt(self, list lids, list pids, double time):
		cdef int i,j
		cdef np.ndarray[double, ndim=1] rel
		cdef LOAD load
		cdef PRESSURE pressure

		# External loads
		for j in lids:
			load = self.mc.loads[j]
			rel = load.Get(time)
			for i in range(load.datRID[0],load.datRID[1]):
				self.datR[i] = rel[i-load.datRID[0]]

		# Pressures
		for j in pids:
			pressure = self.mc.pressures[j]
			rel = pressure.BuildVector(self.dc)
			for i in range(pressure.datRID[0],pressure.datRID[1]):
				self.datR[i] = rel[i-pressure.datRID[0]]

	cpdef ResidualDataCont(self, list cids, int it):
		cdef int i,j
		cdef np.ndarray[double, ndim=1] rel
		cdef CONTACT contact

		# Contacts
		for j in cids:
			contact = self.mc.contacts[j]
			#if it<4:
			contact.FindActiveSlaves(self.dc)
			rel = contact.BuildVector(self.dc)
			for i in range(contact.datRID[0],contact.datRID[1]):
				self.datR[i] = rel[i-contact.datRID[0]]

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
		cdef CONTACT contact
		cdef PRESSURE pressure
		cdef list rowK=[], colK=[], rowR = [], 
		cdef int idK0 = 0, idK1 = 0, idR0=0, idR1=0
		cdef int i, j

		for el in self.mc.elements:
			dof = el.DofID().tolist();

			# Stiffness matrix
			rowK.extend([i for j in range(len(dof)) for i in dof])
			colK.extend([i for i in dof for j in range(len(dof))])
			idK1 = idK0+len(dof)*len(dof);
			el.datKID = [idK0, idK1]
			idK0 = idK1;

			# Residual vector - internal load
			rowR.extend(dof)
			idR1 = idR0+len(dof);
			el.datRID = [idR0, idR1]
			idR0 = idR1;

		for load in self.mc.loads:
			# Residual vector - external load
			dof = load.node.dofID;
			rowR.extend(dof)
			idR1 = idR0+len(dof);
			load.datRID = [idR0, idR1]
			idR0 = idR1;

		for contact in self.mc.contacts:
			# Residual vector - slave and master nodes of contact
			dof = [n.dofID[i] for n in contact.slaveNodes for i in range(self.mc.numDispDofPerNode)];
			dof.extend([n.dofID[i] for n in contact.masterNodes for i in range(self.mc.numDispDofPerNode)]);
			print 'Contact dof: ', dof
			rowK.extend([i for j in range(len(dof)) for i in dof])
			colK.extend([i for i in dof for j in range(len(dof))])
			rowR.extend(dof)
			idR1 = idR0+len(dof)
			idK1 = idK0+len(dof)**2;
			
			# Assigne KID and RID to contact and update idK0 and idR0
			contact.datKID = [idK0, idK1]
			contact.datRID = [idR0, idR1]
			idK0 = idK1;
			idR0 = idR1;

		for pressure in self.mc.pressures:
			# Residual vector - slave and master nodes of contact
			dof = [n.dofID[i] for n in pressure.nodes for i in range(self.mc.numDispDofPerNode)];
			print 'Pressure dof: ', dof
			rowK.extend([i for j in range(len(dof)) for i in dof])
			colK.extend([i for i in dof for j in range(len(dof))])
			rowR.extend(dof)
			idR1 = idR0+len(dof)
			idK1 = idK0+len(dof)**2;
			
			# Assigne KID and RID to contact and update idK0 and idR0
			pressure.datKID = [idK0, idK1]
			pressure.datRID = [idR0, idR1]
			idK0 = idK1;
			idR0 = idR1;

		colR = np.zeros(len(rowR)).tolist()

		return (rowK, colK, rowR, colR)
