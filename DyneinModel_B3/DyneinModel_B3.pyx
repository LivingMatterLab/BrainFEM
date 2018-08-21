from ModelContainer cimport *
from MAT_LE cimport * 
from MAT_VISC cimport * 
from PBAR cimport *
from NodeX cimport *
from CBARX cimport *
from SPC cimport *
from MPC cimport *
from LOAD cimport *
from MT cimport *
from MECH_EL03 cimport *
from Solver cimport *

cimport OutputHelper as oh
cimport ElementHelper as eh

from math import *
import numpy as np
cimport numpy as np

cdef class DyneinModel_B3(ModelContainer):
	def __init__(self):
		super().__init__()
		print "Initializing DyneinModel_B3"

	cpdef BuildModel(self,object p):
		self.numDofPerNode = 3
		self.lAxon = p.lAxon
		self.optionDynein = p.optionDynein

		# Create amplitude
		self.amplitudes.append(Amplitude(np.array([0.,p.tEnd]),np.array([0.,1.])));


		# Compute geometry, which contains
		# xMT 			list of lists. Outer lists contains all MT-lines, inner lists contain [x0,x1,nConn] of each MT element in that line
		# yMT   		list of y coord of each MT-line
		# zMT   		list of z coord of each MT-line
		# crossLinks    list of [MT1_id, MT2_id, x0, x1] of each crossLink
		geom = self.getGeometry(p)
		xMT = geom['xMT']; yMT = geom['yMT']; zMT = geom['zMT']; crossLinks = geom['crossLinks']

		# Create a NEOH and VISC material
		cdef dict matProp = {'E':p.E_MT,'nu':0.4}
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 10;

		matProp['E'] = p.E_Dyn
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 11;

		matProp = {'eta':p.eta_Med}
		self.materials.append(MAT_VISC(matProp))
		self.materials[-1].localID  = 12;

		for m in self.materials:
			print m

		# Create a PBAR properties
		cdef dict propProp = {'area':p.area_MT,'force':0.}
		self.properties.append(PBAR(self.materials[0],propProp))
		self.properties[-1].localID = 100;

		propProp['area'] = p.area_Dyn; propProp['force']= p.force_Dyn
		self.properties.append(PBAR(self.materials[1],propProp))
		self.properties[01].localID = 101;

		propProp = {'area':p.area_MT}
		self.properties.append(PBAR(self.materials[2],propProp))
		self.properties[-1].localID = 102;

		for m in self.properties:
			print m

		# Create mechanism
		self.mechanisms.append(MECH_EL03(p.tCont, p.tDest,p.tCrea, p.minInitStretch_Dyn, p.minStretch_Dyn, p.maxInitStretch_Dyn, p.maxStretch_Dyn,p.activeStretch_Dyn))
		
		# Create nodes and MT elements
		self.MTs = [];				# List of list first and last nodeId for each MT, corresponding to xMT
		nCount = 0;						# Node count
		eCount = 0;						# Element count
		for mtlID in range(len(yMT)):
			for mtID in range(len(xMT[mtlID])):
				x0 = xMT[mtlID][mtID][0]		# Start location of MT
				x2 = xMT[mtlID][mtID][1]		# End location of maximum potentially polymerized MT
				if mtID==0 and x0<1e-6:
					x1 = x2 					# First MT in line, does never polymerize or depolymerize
				elif mtID==len(xMT[mtlID])-1:		# Last MT in line, does never polymerize or depolymerize
					x1 = x2;
				else:
					x1 = np.minimum(x2-p.lMTMax+p.lMT0,x2)				# End location active MT
				nelMT = np.ceil((x2-x0)/p.lStep).astype(int)
				nnodMT = nelMT + 1;

				# Compute nodal locations, last element in MT will generally be shorter than the other ones
				xn = np.concatenate([np.linspace(x0,x0+(nelMT-1)*p.lStep,nelMT),x2*np.ones(1)])

				# Check whether xn(-1)==xn(-2)
				if(np.abs(xn[-1]-xn[-2])<1e-6):
					xn = np.delete(xn,[len(xn)-2])
					nelMT -= 1
					nnodMT -= 1

				# Create the nodes and elements
				nidEndActiveMT = 0
				for nid in range(nnodMT):
					self.nodes.append(NodeX([xn[nid],yMT[mtlID],zMT[mtlID]]));
					self.nodes[-1].localID = nid+nCount;
					self.nodes[-1].dofID = range((nid+nCount)*self.numDofPerNode,(nid+nCount+1)*self.numDofPerNode);
					if(nidEndActiveMT==0 and xn[nid]>=x1):
						nidEndActiveMT = nCount+nid
				
				# Add node ids for this MT to nidMT
				try:
					self.MTs[mtlID].append(MT(nCount,nidEndActiveMT,nCount+nnodMT-1))
				except:
					self.MTs.append([MT(nCount,nidEndActiveMT,nCount+nnodMT-1)])
				if mtID>0: 
					self.MTs[mtlID][-1].mtMinus = self.MTs[mtlID][-2]
					self.MTs[mtlID][-2].mtPlus  = self.MTs[mtlID][-1]

				# Set mechanism and element id of first and last element in this MT
				self.MTs[mtlID][-1].e0 = eCount
				self.MTs[mtlID][-1].e2 = eCount+nelMT-1

				for eid in range(nelMT):
					# Compute initial state
					if(eid+nCount+1>nidEndActiveMT and self.MTs[mtlID][-1].e1==-1):
						self.MTs[mtlID][-1].e1 = eid+eCount-1

					self.elements.append(CBARX([self.nodes[eid+nCount], self.nodes[eid+nCount+1]],self.properties[0]));
					self.elements[-1].localID = eid+eCount;
					self.nodes[eid+nCount].elPlus    = self.elements[-1]
					self.nodes[eid+nCount+1].elMinus = self.elements[-1]

					
					# If this is first element in this MT, and this MT is not the first along its line,
					# then assign corresponding elPlus and elMinus to nodes.
					# These elPlus and elMinus 'connect' two MT
					if(eid==0 and not mtID==0):
						self.nodes[eid+nCount].elMinus = self.nodes[self.MTs[mtlID][mtID-1].n1].elMinus
						self.nodes[self.MTs[mtlID][mtID-1].n1].elPlus   = self.elements[-1]
				
				if self.MTs[mtlID][-1].n1 == self.MTs[mtlID][-1].n2:
					self.MTs[mtlID][-1].e1 = self.MTs[mtlID][-1].e2

				# Increment counts
				nCount += nnodMT;
				eCount += nelMT;

		# Create 2 nodes for storage of deleted cross-links
		self.nodes.append(NodeX([-1.,0.,0.]));
		self.nodes[-1].localID = nCount;
		self.nodes[-1].dofID = range(nCount*self.numDofPerNode,(nCount+1)*self.numDofPerNode);
		nCount+=1
		self.nodes.append(NodeX([-1.,1.,0.]));
		self.nodes[-1].localID = nCount;
		self.nodes[-1].dofID = range(nCount*self.numDofPerNode,(nCount+1)*self.numDofPerNode);
		nCount+=1
		self.storageNodes = [self.nodes[-2],self.nodes[-1]]

		# Create nodes for viscous elements
		nViscNodes = [];				# List of nodeIDs for viscous nodes
		for mtlID in range(len(yMT)):
			self.nodes.append(NodeX([-1.,yMT[mtlID],zMT[mtlID]]));
			self.nodes[-1].localID = nCount;
			self.nodes[-1].dofID = range(nCount*self.numDofPerNode,(nCount+1)*self.numDofPerNode);
			nViscNodes.append(nCount)
			nCount+=1
		
		# Create viscous elements connecting each MT to the clamped side
		initState = NoState
		initTimeToNextEvent = float("inf")
		for mtlID in range(len(yMT)):
			for mtID in range(len(xMT[mtlID])):
				# Create viscous bar element
				self.elements.append(CBARX([self.nodes[nViscNodes[mtlID]], self.nodes[self.MTs[mtlID][mtID].n0]], \
					                       self.properties[2],initState,initTimeToNextEvent));
				self.elements[-1].localID = eCount;
				eCount+=1;
		
		# Create dynein crossLink elements
		self.eidCL0 = eCount
		for csid in range(len(crossLinks)):
			mtl0 = crossLinks[csid][0]			# MTL of first connection point
			mt0  = crossLinks[csid][1]			# MT in this MTL of first connection point
			mtl1 = crossLinks[csid][2]			# MTL of second connection point
			mt1  = crossLinks[csid][3]			# MT in this MTL of second connection point
			x0   = crossLinks[csid][4]			# x of first connection point
			x1   = crossLinks[csid][5]			# x of second connection point

			xmt0 = xMT[mtl0][mt0][0]			# x of start point of first connection MT
			xmt1 = xMT[mtl1][mt1][0]			# x of start point of second connection MT

			
			# Compute node ids of connection points
			nid0 = self.MTs[mtl0][mt0].n0 + np.round((x0-xmt0)/p.lStep).astype(int)
			nid1 = self.MTs[mtl1][mt1].n0 + np.round((x1-xmt1)/p.lStep).astype(int)

			self.elements.append(CBARX([self.nodes[nid0], self.nodes[nid1]],self.properties[1]));	
			self.elements[-1].localID = eCount;
			self.elements[-1].mechanism = self.mechanisms[0];
			self.elements[-1].mechanism.Initialize(self.elements[-1],self);
			eCount+=1;
		self.eidCL1 = eCount-1

		# Compute nodeIDs of nodes at left (x=0) and right (x=lAxon) ends
		nidLeft = [];
		nidRight = [];

		for i in range(len(yMT)):
			nidL = self.MTs[i][0].n0
			if self.nodes[nidL].x<1e-6:
				nidLeft.append(nidL)

			nidR = self.MTs[i][-1].n1;
			if self.nodes[nidR].x>(p.lAxon-1e-6):
				nidRight.append(nidR)
		# Extract node that is loaded
		self.loadNode = self.nodes[nidRight[0]]
		# MPC
		# Constrain right nodes to move same amount in x direction
		for i in range(1,len(nidRight)):
			self.mpc.append(MPC([self.nodes[nidRight[i]],self.loadNode],[0,0],[-1.,1.]))

		# SPC
		# Clamp left nodes in x direction
		for i in range(len(nidLeft)):
			self.spc.append(SPC(self.nodes[nidLeft[i]],[0],0.,self.amplitudes[0]))

		# Clamp storage nodes in x direction
		for n in self.storageNodes:
			self.spc.append(SPC(n,[0],0.,self.amplitudes[0]))

		# Clamp nodes attached to viscous elements in x direction
		for i in nViscNodes:
			self.spc.append(SPC(self.nodes[i],[0],0.,self.amplitudes[0]))

		if p.optionYZ==0:
			# Clamp all nodes in y,z direction
			for n in self.nodes:
				self.spc.append(SPC(n,[1,2],0.,self.amplitudes[0]))
		else:
			# Constrain MTs to have only translation in y-z direction
			for i in range(len(self.MTs)):
				for j in range(len(self.MTs[i])):
					if j==0:
						for s in range(self.MTs[i][j].n0+1,self.MTs[i][j].n1+1):
							self.mpc.append(MPC([self.nodes[s],self.nodes[self.MTs[i][j].n0]],[1,1],[-1.,1.]))
							self.mpc.append(MPC([self.nodes[s],self.nodes[self.MTs[i][j].n0]],[2,2],[-1.,1.]))
					else:
						for s in range(self.MTs[i][j].n0,self.MTs[i][j].n1):
							self.mpc.append(MPC([self.nodes[s],self.nodes[self.MTs[i][j].n1]],[1,1],[-1.,1.]))
							self.mpc.append(MPC([self.nodes[s],self.nodes[self.MTs[i][j].n1]],[2,2],[-1.,1.]))
			# Clamp left and right nodes in y,z direction
			for i in range(len(nidLeft)):
				self.spc.append(SPC(self.nodes[nidLeft[i]],[1,2],0.,self.amplitudes[0]))
			for i in range(len(nidRight)):
				self.spc.append(SPC(self.nodes[nidRight[i]],[1,2],0.,self.amplitudes[0]))

			# Clamp storage nodes in y,z direction
			for n in self.storageNodes:
				self.spc.append(SPC(n,[1,2],0.,self.amplitudes[0]))

			# Clamp nodes attached to viscous elements in y,z direction
			for i in nViscNodes:
				self.spc.append(SPC(self.nodes[i],[1,2],0.,self.amplitudes[0]))

		# Load master node at right in x direction
		self.loads.append(LOAD(self.loadNode,np.array([p.loadExt,0.,0.]),self.amplitudes[0]))
		
		# length of all matrices
		self.SetListLengths()
				
	cpdef getGeometry(self,object p):

		import matplotlib as mpl
		import matplotlib.pyplot as plt
		from mpl_toolkits.mplot3d import Axes3D

		cdef np.ndarray[double,ndim=1]  yMT, zMT
		cdef list xMT, connCS, crossLinks

		# Compute y,z of MicroTubules

		yMT = np.array([0.,p.rInner])
		zMT = np.array([0.,0.])

		# Compute possible cross-links
		#connCS = [[0,1,p.rInner],[1,0,p.rInner]]
		connCS = [[0,1,p.rInner],]	
		nCS = len(connCS)							# Number of potential crosslinks

		# Plot cross section
		plt.figure()
		for i in connCS:
			plt.plot([yMT[i[0]],yMT[i[1]]],[zMT[i[0]],zMT[i[1]]],'k',linewidth=5.0)
		plt.plot(yMT,zMT,'ro',markersize=25)
		plt.plot(p.rAxon*np.cos(np.linspace(0, 2*np.pi, 100)), p.rAxon*np.sin(np.linspace(0, 2*np.pi, 100)),linewidth=5.0)
		plt.axis('equal')
		plt.axis('off')
		plt.xlabel('y')
		plt.ylabel('z')
		plt.savefig(self.folderInput+'/cs.png',format='png')

		# Compute length of MT at x=0
		xMT = [[[0.,p.lMT0,0],],[[p.lAxon-p.lMT0,p.lAxon,0],]]

		# Compute cross links
		crossLinks = [] # Every entry contains MTL1, MT1, MTL2, MT2, x0, x1
		for i in range(1,int(p.lAxon/p.dlLink)):
			x0 = i*p.dlLink
			csID = np.random.randint(0,nCS)
			x1 = x0-connCS[csID][2]*np.tan(p.thLink)

			# Check that there is no gap at either one of the MT
			xMT0 = np.asarray(xMT[connCS[csID][0]])
			xMT1 = np.asarray(xMT[connCS[csID][1]])

			try:
				id01 = [ n for n,i in enumerate(xMT0[:,0]) if i<=x0][-1]
				id02 = [ n for n,i in enumerate(xMT0[:,1]) if i>=x0][0]

				id11 = [ n for n,i in enumerate(xMT1[:,0]) if i<=x1][-1]
				id12 = [ n for n,i in enumerate(xMT1[:,1]) if i>=x1][0]

				# if id01==id02 and id11==id12, both endpoints of the crosslink are attached
				if(id01==id02 and id11==id12):
					crossLinks.append([connCS[csID][0],id01,connCS[csID][1],id11,x0,x1])

					# Increment counter for cross links on these MT
					xMT[connCS[csID][0]][id01][2] +=1
					xMT[connCS[csID][1]][id11][2] +=1

			except:
				pass


		# Plot all MT and crossLinks in 3D
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		colorMT = 'bgrcmy'
		# MT
		for i in range(2):
			for j in range(len(xMT[i])):
				ax.plot([ xMT[i][j][0],xMT[i][j][1]],[yMT[i],yMT[i]],[zMT[i],zMT[i]],color=colorMT[j%len(colorMT)])
		# CrossLinks
		for i in range(len(crossLinks)):
			ax.plot([ crossLinks[i][4]     ,crossLinks[i][5] ],\
					[ yMT[crossLinks[i][0]],yMT[crossLinks[i][2]] ],\
					[ zMT[crossLinks[i][0]],zMT[crossLinks[i][2]] ],'k',linewidth=0.5)

		for direction in (-1, 1):
			for point in np.diag(direction * p.lAxon * np.array([0,1,1])):
				ax.plot([point[0]], [point[1]], [point[2]], 'w')
		plt.savefig(self.folderInput+'/axon1.eps',format='eps')
		ax.view_init(elev=0, azim=0)
		plt.savefig(self.folderInput+'/axon2.eps',format='eps')
		ax.view_init(elev=0, azim=-90)
		plt.savefig(self.folderInput+'/axon3.eps',format='eps')

		# Compute the number of cross-links for each MT element
		nCS_El = []
		nCS_0  = 0;		# Number of elements with 0 cross-links
		for i in xMT:
			for j in i:
				nCS_El.append(j[2])
				if(j[2]==0):
					nCS_0+=1
		print 'The axon model has '+ str(nCS_0) + ' unconnected elements'

		fig = plt.figure()
		plt.hist(nCS_El, 50, facecolor='green', alpha=0.75)
		plt.title((str(nCS_0)+' unconnected elements'))
		plt.xlabel('# of cross links at MT')
		plt.ylabel('frequency')
		plt.savefig(self.folderInput+'/cs_hist.png',format='png')

		return {'xMT':xMT, 'yMT':yMT,'zMT':zMT,'crossLinks':crossLinks}

	cpdef WriteStepOutput(self, Solver s):
		import params as p
		cdef str strWrite
		cdef Node n0, n1;
		cdef int i, j
		cdef list n0Dof, n1Dof

		# MT.txt contains location of all MT in this step
		strWrite = "Time " + '{: <15.6f}'.format(s.dc.time) +"\n"
		strWrite+= '{: <15}'.format('MT') + '{: <15}'.format('N1.X') + '{: <15}'.format('N1.Y') + '{: <15}'.format('N1.Z')+ '{: <15}'.format('N2.X') + '{: <15}'.format('N2.Y') + '{: <15}'.format('N2.Z') +'\n'

		MTcount = 1
		for i in range(len(self.MTs)):
			for mt in self.MTs[i]:
				n0 = self.nodes[mt.n0]
				n1 = self.nodes[mt.n1]


				n0Dof = n0.Dof(s.dc)
				n1Dof = n1.Dof(s.dc)
				while len(n0Dof)<len(n0.loc): n0Dof.append(0.0);
				while len(n1Dof)<len(n1.loc): n1Dof.append(0.0);
				
				strWrite+= '{: <15d}'.format(MTcount) + \
				           '{: <15.6f}'.format(n0.x+n0Dof[0]) + \
				           '{: <15.6f}'.format(n0.y+n0Dof[1]) + \
				           '{: <15.6f}'.format(n0.z+n0Dof[2]) + \
				           '{: <15.6f}'.format(n1.x+n1Dof[0]) + \
				           '{: <15.6f}'.format(n1.y+n1Dof[1]) + \
				           '{: <15.6f}'.format(n1.z+n1Dof[2]) + \
				           '\n'

				MTcount +=1
		oh.WriteToOutput(self,'MT.txt',strWrite)


		# CL.txt contains location of all crosslinks in this step
		strWrite = "Time " + '{: <15.6f}'.format(s.dc.time) +"\n"
		strWrite+= '{: <15}'.format('Element') + '{: <15}'.format('N1.X') + '{: <15}'.format('N1.Y') + '{: <15}'.format('N1.Z')+ '{: <15}'.format('N2.X') + '{: <15}'.format('N2.Y') + '{: <15}'.format('N2.Z') +'\n'

		CLcount = 0
		for i in range(self.eidCL0,self.eidCL1+1):
			n0 = self.elements[i].nodes[0]
			n1 = self.elements[i].nodes[1]

			# Check if this cross link is released
			if n0.localID == self.storageNodes[0].localID: continue
			
			n0Dof = n0.Dof(s.dc)
			n1Dof = n1.Dof(s.dc)
			while len(n0Dof)<len(n0.loc): n0Dof.append(0.0);
			while len(n1Dof)<len(n1.loc): n1Dof.append(0.0);
			
			strWrite+= '{: <15d}'.format(self.elements[i].localID) + \
			           '{: <15.6f}'.format(n0.x+n0Dof[0]) + \
			           '{: <15.6f}'.format(n0.y+n0Dof[1]) + \
			           '{: <15.6f}'.format(n0.z+n0Dof[2]) + \
			           '{: <15.6f}'.format(n1.x+n1Dof[0]) + \
			           '{: <15.6f}'.format(n1.y+n1Dof[1]) + \
			           '{: <15.6f}'.format(n1.z+n1Dof[2]) + \
			           '\n'

			CLcount +=1
		oh.WriteToOutput(self,'CL.txt',strWrite)

		# FD.txt contains force-displacement at axon tip
		if s.dc.step==0:
			strWrite = '{: <16s}'.format('Step') + '{: <16s}'.format('Time') + \
					   '{: <16s}'.format('Displ.') +  '{: <16s}'.format('Stretch') + \
					   '{: <16s}'.format('Force') +  '{: <16s}'.format('Piola1')  + \
					   '{: <16s}'.format('# of CL') + '\n'
		else:
			strWrite = '';


		strWrite+= '{: <16d}'.format(s.dc.step) + \
		           '{: <16.6f}'.format(s.dc.time) + \
		           '{: <16.6e}'.format(self.loadNode.Dof(s.dc)[0]) + \
		           '{: <16.6f}'.format(1.+self.loadNode.Dof(s.dc)[0]/p.lAxon) + \
		           '{: <16.6e}'.format(min(1.,s.dc.time/p.tLoad)*p.loadExt) + \
		           '{: <16.6e}'.format(min(1.,s.dc.time/p.tLoad)*p.loadExt/pi/p.rAxon**2) + \
		           '{: <16d}'.format(CLcount)

		oh.WriteToOutput(self,'FD.txt',strWrite)

	cpdef UpdateModel(self, Solver s):
		cdef Element el
		cdef MT mt
		cdef int i
		cdef int nDestroy, nCreate, nC
		cdef double dist0, dist,frac, randMT
		cdef np.ndarray[np.int_t, ndim=1] dof

		# Clear list of restores
		self.elRestore = [];
		self.mtRestore = [];

		# Decrement timeToNextEvent in all elements and MT
		for el in self.elements:
			el.timeToNextEvent-=s.dc.dt

		for i in range(len(self.MTs)):
			for mt in self.MTs[i]:
				mt.timeToNextEvent-=s.dc.dt

		
		# First check whether Microtubule event has to happen
		for i in range(len(self.MTs)):
			for mt in self.MTs[i]:
				try:
					mt.mechanism.Apply(mt,self,s)
				except:
					pass

		# Check for events in cross-link elements
		nDestroy = 0;
		nCreate  = 0;
		for el in self.elements:

			# Apply mechanism
			if not el.mechanism:
				nC = 0;
			else:
				nC = el.mechanism.Apply(el,self,s)
			
			if(nC>0): 
				nCreate+=nC; 
			else: 
				nDestroy-=nC;

		oh.WriteToLog(self,str(nDestroy)+" crosslinks destroyed.")
		oh.WriteToLog(self,str(nCreate) +" crosslinks created.")


	cpdef RestoreModel(self, Solver s):
		cdef Element el
		cdef int ndof, countK, rid,i, j
		cdef np.ndarray[np.int_t, ndim=1] dof
		cdef dict res		

		for res in reversed(self.elRestore):
			el       			= res['el']
			el.state 			= res['state']
			el.timeToNextEvent 	= res['timeToNextEvent'] 

			try:
				el.restLength       = res['restLength']
				el.nodes 			= [self.nodes[i] for i in res['nodes']]
				el.dummyNodes 		= [self.nodes[i] for i in res['dummyNodes']]


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
			except:
				pass

		for res in reversed(self.mtRestore):
			mt                 = res['mt']
			mt.n1              = res['n1']
			mt.e1              = res['e1']
			mt.state           = res['state']
			mt.timeToNextEvent = res['timeToNextEvent']

			# Recompute all element pointers related to this MT
			for eid in range(mt.e0,mt.e2+1):
				el = self.elements[eid]
				el.nodes[0].elPlus = el
				el.nodes[1].elMinus = el
			self.elements[mt.e1].nodes[1].elPlus = self.elements[mt.mtPlus.e0];
			self.nodes[mt.mtPlus.n0].elMinus = self.elements[mt.e1]

		# Clear restore
		self.elRestore = [];
		self.mtRestore = [];

		# Finally restore timeToNextEvent in all elements and MT
		for el in self.elements:
			el.timeToNextEvent+=s.dc.dt
		for i in range(len(self.MTs)):
			for mt in self.MTs[i]:
				mt.timeToNextEvent+=s.dc.dt

		oh.WriteToLog(self,"Creation and destruction of crosslinks restored.\n")
