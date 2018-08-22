from MTwithCortex_B3 cimport *
from Amplitude cimport *
from MAT_LE cimport * 
from MAT_VISC cimport * 
from PBAR cimport *
from PBAR1 cimport *
from NodeX cimport *
from CBARX2 cimport *
from SPC cimport *
from MPC cimport *
from LOAD cimport *
from MT cimport *
from MECH_EL01 cimport *
from MECH_EL02 cimport *
from MECH_EL03 cimport *
from MECH_EL04 cimport *
from MECH_MT02 cimport *
from Solver cimport *

cimport OutputHelper as oh
cimport ElementHelper as eh

from math import *
import numpy as np
cimport numpy as np

cdef class MTwithCortex2_B3(MTwithCortex_B3):
	def __init__(self):
		super().__init__()
		print "Initializing MTwithCortex2_B3"

	cpdef BuildModel(self,object p):
		self.numDofPerNode = 1
		self.lAxon = p.lAxon
		self.optionDynein = p.optionDynein
		
		#  ================================================  #
		#  ---------------  G E O M E T R Y  --------------  #
		#  ------------------------------------------------  #
		
		# Compute geometry, which contains
		# xMT 			list of lists. Outer lists contains all MT-lines, inner lists contain [x0,x1,nConn] of each MT element in that line
		# yMT   		list of y coord of each MT-line
		# zMT   		list of z coord of each MT-line
		# crossLinks    list of [MT1_id, MT2_id, x0, x1] of each crossLink
		geomMT = self.getGeometry(p)
		
		# Compute cortex geometry, which contains
		# xAc 			list of lists. Outer lists contains all Actin-lines, inner lists contain [x0,x1,nConn] of each Actin element in that line
		# yAc   		list of y coord of each Actin-line
		# zAc   		list of z coord of each Actin-line
		# crossLinks    list of [AcL0, Ac0, AcL1, Ac1, x0, x1] OR [AcL0, Ac0, -1, Ring1, x0, th1] of each crosslink
		# xAcRing		array of x-coordinates of each actin filament ring
		geomAc = self.getCortexGeometry(p)
		#geomAc = self.getSimpleCortexGeometry(p)

		# Compute crosslinks between actin and MT, which contains:
		# crossLinks    list of [MTL0, MT0, AcL1, Ac1, x0, x1] of each crosslink
		geomCL_Ac_Mt = self.getCL_MT_Ac(p,geomMT,geomAc)
		#geomCL_Ac_Mt = {'crossLinks':[]}

		#  ================================================  #
		#    M A T E R I A L S  A N D  P R O P E R T I E S   #
		#  ------------------------------------------------  #
		# Create amplitude
		#self.amplitudes.append(Amplitude(np.array([0.,p.tLoad,p.tEnd]),np.array([0.,1.,1.])));
		self.amplitudes.append(Amplitude(np.array([0.,p.tLoad,p.tLoad+p.tHold,p.tLoad+p.tHold+p.tFinalLoad]),np.array([0.,1.,1.,2.])));



		# Create a NEOH and VISC material
		cdef dict matProp = {'E':p.E_MT,'nu':0.4}
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 10;				# MT - elastic

		matProp['E'] = p.E_Dyn
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 11;				# Dynein - elastic

		matProp['E'] = p.E_Ac
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 12;				# Actin - elastic
		
		matProp['E'] = p.E_Myo
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 13;				# Myosin - elastic

		matProp['E'] = p.E_MT_Ac
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 14;				# MT Actin cl - elastic

		matProp['E'] = p.E_GC
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 15;				# Growth cone conn. - elastic

		matProp['E'] = p.E_Med
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 16;				# Medium - elastic

		matProp = {'eta':p.eta_Ac}					
		self.materials.append(MAT_VISC(matProp))
		self.materials[-1].localID  = 17;				# Actin - viscous

		matProp = {'eta':p.eta_Med}
		self.materials.append(MAT_VISC(matProp))
		self.materials[-1].localID  = 18;				# Medium - viscous

		for m in self.materials:
			print m

		# Create a PBAR properties
		cdef dict propProp = {'area':p.area_MT,'force':0.}
		self.properties.append(PBAR(self.materials[0],propProp))
		self.properties[-1].localID = 100;				# MT - elastic

		propProp['area'] = p.area_Dyn; propProp['force']= p.force_Dyn
		self.properties.append(PBAR1(self.materials[1],propProp,self.amplitudes[0]))
		self.properties[-1].localID = 101;				# Dynein - elastic

		propProp['area'] = p.area_Ac; propProp['force']= p.force_Ac
		self.properties.append(PBAR1(self.materials[2],propProp,self.amplitudes[0]))
		self.properties[-1].localID = 202;				# Actin - elastic

		propProp['area'] = p.area_Myo; propProp['force']= p.force_Myo
		self.properties.append(PBAR1(self.materials[3],propProp,self.amplitudes[0]))
		self.properties[-1].localID = 203;				# Myosin - elastic

		propProp['area'] = p.area_MT_Ac; propProp['force']= 0.
		self.properties.append(PBAR1(self.materials[4],propProp,self.amplitudes[0]))
		self.properties[-1].localID = 304;				# MT_Ac - elastic

		propProp['area'] = p.area_GC; propProp['force']= 0.
		self.properties.append(PBAR1(self.materials[5],propProp,self.amplitudes[0]))
		self.properties[-1].localID = 105;				# MT_Ac - elastic

		propProp = {'area':p.area_MT}
		self.properties.append(PBAR(self.materials[6],propProp))
		self.properties[-1].localID = 406; 				# Medium - elastic

		propProp['area'] = p.area_Ac;
		self.properties.append(PBAR(self.materials[7],propProp))
		self.properties[-1].localID = 207;				# Actin - viscous

		propProp = {'area':p.area_MT}
		self.properties.append(PBAR(self.materials[8],propProp))
		self.properties[-1].localID = 408;				# Medium - viscous

		for m in self.properties:
			print m

		# Create mechanism
		#self.mechanisms.append(MECH_EL02(p.tCrea, p.tDest, p.maxInitStretch_Dyn, p.maxStretch_Dyn))
		self.mechanisms.append(MECH_EL03(p.tCont_Dyn, p.tDest_Dyn,p.tCrea_Dyn,\
		                                 p.minInitStretch_Dyn, p.minStretch_Dyn,\
		                                 p.maxInitStretch_Dyn, p.maxStretch_Dyn, \
		                                 p.activeStretch_Dyn))
		self.mechanisms.append(MECH_EL02(p.tCrea_Tau,p.tDest_Tau,\
		                                 p.maxInitStretch_Dyn, p.maxStretch_Dyn))
		self.mechanisms.append(MECH_MT02(p.tMTpoly,p.tMTstat,p.tMTdepoly,\
			                             p.MTpolyRate,p.MTdepolyRate, \
			                             p.dlGC/2.,p.stabDistalMT,
			                             p.fracLMT,p.lMT0/p.lMTMax*19.))
		self.mechanisms.append(MECH_EL04(p.tCont_Myo, p.tDest_Myo,p.tCrea_Myo,\
		                                 p.minInitStretch_Myo, p.minStretch_Myo,\
		                                 p.maxInitStretch_Myo, p.maxStretch_Myo,\
		                                 p.activeStretch_Myo))
		self.mechanisms.append(MECH_EL04(p.tCont_Myo, p.tDest_Myo,p.tCrea_Myo,\
		                                 p.minInitStretch_Myo, p.minStretch_Myo,\
		                                 p.maxInitStretch_Myo, p.maxStretch_Myo,\
		                                 1.0))
		self.mechanisms.append(MECH_EL02(p.tCrea_MT_Ac, p.tDest_MT_Ac, \
			                             p.maxInitStretch_Dyn, p.maxStretch_Dyn))
		

		#  ================================================  #
		#  ------  N O D E S  A N D  E L E M E N T S  -----  #
		#  ------------------------------------------------  #
		self.MTs = [];				# List of list first and last nodeId for each MT, corresponding to xMT
		nCount = 0;						# Node count
		eCount = 0;						# Element count

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


		# Create nodes and MT elements
		if p.includeMT:
			nCount, eCount, nMTViscNodes = self.CreateMT(p,geomMT, nCount, eCount)

		# Create nodes and elements for longitudinal actin filaments
		self.Acs = [];					# List of list first and last nodeId for each Actin, corresponding to xAc
		if p.includeAc:
			nCount, eCount, nAcViscNode, nidAcRing0, nidAcRing1 = self.CreateActin(p,geomAc, nCount, eCount)

		"""
		for acl in self.Acs:
			for ac in acl:
				for eid in range(ac.e0,ac.e2):
					print self.elements[eid]
				print '\n'
			print '\n'
		"""

		# Create dynein crossLink elements among MT
		if p.includeMT:
			nCount, eCount = self.CreateMTCrosslinks(p,geomMT,nCount,eCount);

		# Create myosin crosslink elements within the cortex
		if p.includeAc:
			nCount, eCount = self.CreateActinCrosslinks(p,geomAc,nidAcRing0, nidAcRing1, nCount,eCount);

		# Create coupling crossLink elements between MT and cortex
		if p.includeMT and p.includeAc and p.includeAcMT:
			nCount, eCount = self.CreateAcMTCrosslinks(p,geomMT,geomAc,geomCL_Ac_Mt,nCount,eCount);

		# Create connections elements between MT and growth cone
		#nCount, eCount = self.CreateGCConnections(p,geomMT,nCount,eCount);

		#  ================================================  #
		#  ----  B O U N D A R Y  C O N D I T I O N S  ----  #
		#  ------------------------------------------------  #
		# Compute nodeIDs of nodes at left (x=0) and right (x=lAxon) ends
		nidLeft = [];
		nidRightMT = [];
		nidRightAc = [];
		for n in self.nodes:
			if np.abs(n.x)<1.e-6 and np.abs(n.y)<2*p.rAxon and np.abs(n.z)<2*p.rAxon:
				nidLeft.append(n.localID)
			if np.abs(n.x-p.lAxon)<1.e-6:
				if np.sqrt(np.abs(n.y)**2+np.abs(n.z)**2)<0.99*p.rAxon:
					nidRightMT.append(n.localID)
				elif np.sqrt(np.abs(n.y)**2+np.abs(n.z)**2)<2.*p.rAxon:
					nidRightAc.append(n.localID)
		
		# Extract node that is loaded
		#self.loadNodeMT = self.nodes[nidRightMT[0]]

		# MPC
		# Constrain right nodes to move same amount in x direction
		#for i in range(1,len(nidRightMT)):
		#	self.mpc.append(MPC([self.nodes[nidRightMT[i]],self.loadNodeMT],[0,0],[-1.,1.]))
		if p.includeAc and p.includeMT:
			for i in range(len(nidRightAc)):
				self.mpc.append(MPC([self.nodes[nidRightAc[i]],self.growthConeNode],[0,0],[-1.,1.]))
		elif p.includeAc:
			self.growthConeNode = self.nodes[nidRightAc[0]]
			for i in range(1,len(nidRightAc)):
				self.mpc.append(MPC([self.nodes[nidRightAc[i]],self.growthConeNode],[0,0],[-1.,1.]))

		# Constrain all actin rings to move same amount in x direction
		if p.includeAc:
			for i in range(len(nidAcRing0)):
				for nid in range(nidAcRing0[i],nidAcRing1[i]):
					self.mpc.append(MPC([self.nodes[nid],self.nodes[nidAcRing1[i]]],[0,0],[-1.,1.]))
		
		# SPC
		# Clamp left nodes in x direction
		for i in range(len(nidLeft)):
			self.spc.append(SPC(self.nodes[nidLeft[i]],[0],0.,self.amplitudes[0]))

		# Clamp storage nodes in x direction
		for n in self.storageNodes:
			self.spc.append(SPC(n,[0],0.,self.amplitudes[0]))

		# Clamp viscNodes in x direction
		if p.includeMT:
			for nid in nMTViscNodes:
				self.spc.append(SPC(self.nodes[nid],[0],0.,self.amplitudes[0]))
		if p.includeAc:
			self.spc.append(SPC(self.nodes[nAcViscNode],[0],0.,self.amplitudes[0]))
		
		# Load master node at right in x direction
		if p.includeMT:
			nid = self.growthConeNode.localID;
		else:
			nid = nidRightAc[0]
		if p.optionBC==0:	
			self.loads.append(LOAD(self.nodes[nid],np.array([p.loadExt,0.,0.]),self.amplitudes[0]))
		else:
			self.spc.append(SPC(self.nodes[nid],[0],p.dispExt,self.amplitudes[0]))

		# length of all matrices
		self.SetListLengths()

	cpdef CreateMT(self, object p, dict geom, int nCount, int eCount):
		xMT = geom['xMT']; yMT = geom['yMT']; zMT = geom['zMT'];

		for mtlID in range(len(yMT)):
			for mtID in range(len(xMT[mtlID])):
				x0 = xMT[mtlID][mtID][0]		# Start location of MT
				x2 = xMT[mtlID][mtID][1]		# End location of maximum potentially polymerized MT
				if mtID==0 and x0<1e-6:
					x1 = x2 					# First MT in line, does never polymerize or depolymerize
				elif mtID==len(xMT[mtlID])-1:		# Last MT in line cannot extend beyond lAxon+dlGC/2.
					x1 = np.minimum(x2,np.minimum(x0+p.lMT0,p.lAxon+p.dlGC/2.))
				else:
					x1 = np.minimum(x2-p.lMTMax+p.lMT0,x2)				# End location active MT
				nelMT = np.ceil((x2-x0)/p.dl_MT).astype(int)
				nnodMT = nelMT + 1;

				# Compute nodal locations, last element in MT will generally be shorter than the other ones
				xn = np.concatenate([np.linspace(x0,x0+(nelMT-1)*p.dl_MT,nelMT),x2*np.ones(1)])

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
				self.MTs[mtlID][-1].mechanism = self.mechanisms[2]
				self.MTs[mtlID][-1].e0 = eCount
				self.MTs[mtlID][-1].e2 = eCount+nelMT-1

				for eid in range(nelMT):
					# Compute initial state
					if(eid+nCount+1>nidEndActiveMT and self.MTs[mtlID][-1].e1==-1):
						self.MTs[mtlID][-1].e1 = eid+eCount-1

					self.elements.append(CBARX2([self.nodes[eid+nCount], self.nodes[eid+nCount+1]],self.properties[0]));
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

		# Create node for growth cone
		self.nodes.append(NodeX([p.lAxon+p.dlGC,0.,0.]));
		self.nodes[-1].localID = nCount;
		self.nodes[-1].dofID = range(nCount*self.numDofPerNode,(nCount+1)*self.numDofPerNode);
		nCount+=1
		self.growthConeNode = self.nodes[-1]
		
		# Create connections elements between MT and growth cone
		for mtl in self.MTs:
			self.elements.append(CBARX2([self.nodes[mtl[-1].n1], self.growthConeNode],self.properties[5]));	
			self.elements[-1].localID = eCount;
			mtl[-1].eGC = self.elements[-1].localID
			eCount+=1;

		# Initialize mt mechanisms
		for mtlID in range(len(yMT)):
			for mtID in range(len(xMT[mtlID])):
				self.MTs[mtlID][mtID].mechanism.Initialize(self.MTs[mtlID][mtID],self)

		# Create nodes for viscous elements
		nViscNodes = [];				# List of nodeIDs for viscous nodes
		for mtlID in range(len(yMT)):
			self.nodes.append(NodeX([-1.,yMT[mtlID],zMT[mtlID]]));
			self.nodes[-1].localID = nCount;
			self.nodes[-1].dofID = range(nCount*self.numDofPerNode,(nCount+1)*self.numDofPerNode);
			nViscNodes.append(nCount)
			nCount+=1

		# Compute total MT length and set as goalLMT to mechanicm
		#self.currLMT = self.CalcLMT(None)
		#self.MTs[0][0].mechanism.goalLMT = self.currLMT
		
		# Create viscous elements connecting each MT to the clamped side
		initState = NoState
		initTimeToNextEvent = float("inf")
		for mtlID in range(len(yMT)):
			for mtID in range(len(xMT[mtlID])):
				# Create viscous bar element
				self.elements.append(CBARX2([self.nodes[nViscNodes[mtlID]], self.nodes[self.MTs[mtlID][mtID].n0]], \
					                       self.properties[8],initState,initTimeToNextEvent));
				self.elements[-1].localID = eCount;
				eCount+=1;		

		return (nCount, eCount, nViscNodes)

	cpdef CreateMTCrosslinks(self, object p, dict geom, int nCount, int eCount):
		xMT = geom['xMT']; crossLinks = geom['crossLinks']

		self.eidMTCL0 = eCount
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
			nid0 = self.MTs[mtl0][mt0].n0 + np.round((x0-xmt0)/p.dl_MT).astype(int)
			nid1 = self.MTs[mtl1][mt1].n0 + np.round((x1-xmt1)/p.dl_MT).astype(int)

			self.elements.append(CBARX2([self.nodes[nid0], self.nodes[nid1]],self.properties[1]));	
			self.elements[-1].localID = eCount;
			if np.random.rand()<p.fracForce_Dyn:
				self.elements[-1].mechanism = self.mechanisms[0];
			else:
				self.elements[-1].mechanism = self.mechanisms[1];
			self.elements[-1].mechanism.Initialize(self.elements[-1],self);
			eCount+=1;
		self.eidMTCL1 = eCount-1

		return (nCount, eCount)

	cpdef CreateActin(self, object p, dict geom, int nCount, int eCount):
		xAc = geom['xAc']; yAc = geom['yAc']; zAc = geom['zAc'];
		xAcRing = geom['xAcRing']; 

		for aclID in range(len(yAc)):
			acl = xAc[aclID]
			for acID in range(len(acl)):
				ac = acl[acID]
				x0 = ac[0]		# Start x location of Actin
				x2 = ac[1]		# End x location Actin
				x1 = x2
				nelAc = np.ceil((x2-x0)/p.dl_Ac).astype(int)
				nnodAc = nelAc + 1;

				# Compute dx as distance between two nodes
				dx = (x2-x0)/nelAc
				
				# Create the nodes
				xi = x0
				for nid in range(nnodAc):
					self.nodes.append(NodeX([xi,yAc[aclID],zAc[aclID]]));
					self.nodes[-1].localID = nid+nCount;
					self.nodes[-1].dofID = range((nid+nCount)*self.numDofPerNode,(nid+nCount+1)*self.numDofPerNode);
					xi += dx

				# Add segment as MT (actin has minus sides on both ends, and plus in center)
				try:
					self.Acs[aclID].append(MT(nCount,nCount+nnodAc-1,nCount+nnodAc-1))
				except:
					self.Acs.append([MT(nCount,nCount+nnodAc-1,nCount+nnodAc-1)])
				if acID>0: 
					self.Acs[aclID][-1].mtMinus  = self.Acs[aclID][-2]
					self.Acs[aclID][-2].mtMinus  = self.Acs[aclID][-1]
				self.Acs[aclID][-1].e0 = eCount
				self.Acs[aclID][-1].e2 = eCount+nelAc-1
				self.Acs[aclID][-1].e1 = self.Acs[aclID][-1].e2

				# Create elements
				for eid in range(nelAc):
					self.elements.append(CBARX2([self.nodes[eid+nCount], self.nodes[eid+nCount+1]],self.properties[7]));
					self.elements[-1].state = Actin;
					self.elements[-1].localID = eid+eCount;
					if eid<=nelAc/2:
						self.nodes[eid+nCount].elPlus    = self.elements[-1]
						self.nodes[eid+nCount+1].elMinus = self.elements[-1]
					else:
						self.nodes[eid+nCount].elMinus    = self.elements[-1]
						self.nodes[eid+nCount+1].elPlus = self.elements[-1]
					if eid==nelAc/2+1:
						self.nodes[eid+nCount].elMinus = self.elements[-2]
						self.nodes[eid+nCount].elPlus = self.elements[-1]

					# Update polarity of elements
					if self.elements[-1].nodes[0].elPlus == self.elements[-1]:
						self.elements[-1].minusID = 0;
						self.elements[-1].plusID = 1;
					else:
						self.elements[-1].minusID = 1;
						self.elements[-1].plusID = 0;

					# If this is first element in this Actin filament, 
					# and this Actin filament is not the first along its line,
					# then assign corresponding elPlus and elMinus to nodes.
					# These elPlus and elMinus 'connect' two Actin filaments
					if(eid==0 and not acID==0):
						self.nodes[eid+nCount].elMinus = self.nodes[self.Acs[aclID][acID-1].n1].elPlus
						self.nodes[self.Acs[aclID][acID-1].n1].elMinus   = self.elements[-1]

				# Increment counts
				nCount += nnodAc;
				eCount += nelAc;	

		# Create nodes and elements for circular actin filaments
		nR = np.ceil(2*pi*p.rAxon/p.dl_Ac).astype(int)
		thR = np.linspace(0,2*pi,nR+1); thR = np.delete(thR,-1)
		yR = p.rAxon*p.rFracActin*np.cos(thR); 
		zR = p.rAxon*p.rFracActin*np.sin(thR)

		nidRing0 = []; nidRing1 = []
		for i in range(len(xAcRing)):
			if i%2==0:
				yR/=p.rFracActin;
				zR/=p.rFracActin;
			else:
				yR*=p.rFracActin;
				zR*=p.rFracActin;

			# Create nodes
			for nid in range(nR):
				self.nodes.append(NodeX([xAcRing[i],yR[nid],zR[nid]]));
				self.nodes[-1].localID = nid+nCount;
				self.nodes[-1].dofID = range((nid+nCount)*self.numDofPerNode,(nid+nCount+1)*self.numDofPerNode);

			# Create elements
			for eid in range(nR):
				self.elements.append(CBARX2([self.nodes[eid+nCount], self.nodes[(eid+1)%nR+nCount]],self.properties[2]));
				self.elements[-1].state = Actin;
				self.elements[-1].localID = eid+eCount;
			
			# Add element ids for this ring
			nidRing0.append(nCount)
			nidRing1.append(nCount+nR-1)

			# Increment counts
			nCount += nR;
			eCount += nR;

		# Create node far away to attach visc elements representing the medium
		self.nodes.append(NodeX([-10*p.lAxon,0.,0.]));
		self.nodes[-1].localID = nCount;
		self.nodes[-1].dofID = range(nCount*self.numDofPerNode,(nCount+1)*self.numDofPerNode);
		nViscNode = nCount
		nCount += 1

		# Attach one node in each actin ring to the viscNodes.
		for nid in nidRing0:
			# Create CBAR elements for supports
			self.elements.append(CBARX2([self.nodes[nid], self.nodes[nViscNode]],self.properties[8]));
			self.elements[-1].localID = eCount;
			eCount+=1

		# Attach end points of each longitudinal actin filaments to the viscNodes.
		for acl in self.Acs:
			for ac in acl:
				self.elements.append(CBARX2([self.nodes[ac.n0], self.nodes[nViscNode]],self.properties[8]));
				self.elements[-1].localID = eCount;
				eCount+=1

		return (nCount, eCount, nViscNode, nidRing0, nidRing1)

	cpdef CreateActinCrosslinks(self,object p,dict geom, list nidRing0, list nidRing1, int nCount, int eCount):
		crossLinks = geom['crossLinks']

		nR = np.ceil(2*pi*p.rAxon/p.dl_Ac).astype(int)
		thR = np.linspace(0,2*pi,nR+1); thR = np.delete(thR,-1)
		
		# Create crossLink elements within actin
		self.eidAcCL0 = eCount
		for cl in crossLinks:

			# First connection point
			acl0 = cl[0]			# AcL of first connection point
			ac0  = cl[1]			# Ac in this AcL of first connection point
			x0   = cl[4]			# x of first connection point

			n0 = self.nodes[self.Acs[acl0][ac0].n0]		# First node at this actin filament
			n1 = self.nodes[self.Acs[acl0][ac0].n1]		# Last node at this actin filament
			nid0 = n0.localID+np.round((n1.localID-n0.localID)*(x0-n0.x)/(n1.x-n0.x)).astype(int)
			
			if cl[2]==-1:
				r1 = cl[3]
				th1 = cl[5]

				n0 = self.nodes[nidRing0[r1]]		# First node at this actin ring
				nid1 = n0.localID+np.round(th1/(thR[1]-thR[0])).astype(int)%nR
			else:
				# Cross link between two longitudinal actin filaments
				acl1 = cl[2]			# AcL of second connection point
				ac1  = cl[3]			# Ac in this AcL of second connection point
				x1   = cl[5]			# x of second connection point

				n0 = self.nodes[self.Acs[acl1][ac1].n0]		# First node at this actin filament
				n1 = self.nodes[self.Acs[acl1][ac1].n1]		# Last node at this actin filament
				nid1 = n0.localID+np.round((n1.localID-n0.localID)*(x1-n0.x)/(n1.x-n0.x)).astype(int)
				

			self.elements.append(CBARX2([self.nodes[nid0], self.nodes[nid1]],self.properties[3]));	
			self.elements[-1].localID = eCount;
			self.elements[-1].state = MotorAttachedTaut;
			if not cl[2]==-1:
				if np.random.rand()<p.fracForce_Myo:
					self.elements[-1].mechanism = self.mechanisms[3];
				else:
					self.elements[-1].mechanism = self.mechanisms[4];
				self.elements[-1].mechanism.Initialize(self.elements[-1],self);
			eCount+=1;
		self.eidAcCL1 = eCount-1

		return (nCount, eCount)

	cpdef CreateAcMTCrosslinks(self, object p, dict geomMT, dict geomAc,\
	                           dict geomCL_Ac_Mt, int nCount, int eCount):
		xMT = geomMT['xMT'];
		xAc = geomAc['xAc']; 
		crossLinks = geomCL_Ac_Mt['crossLinks']

		self.eidMTAcCL0 = eCount
		for csid in range(len(crossLinks)):
			mtl0 = crossLinks[csid][0]			# MTL of first connection point
			mt0  = crossLinks[csid][1]			# MT in this MTL of first connection point
			acl1 = crossLinks[csid][2]			# AcL of second connection point
			ac1  = crossLinks[csid][3]			# Ac in this AcL of second connection point
			x0   = crossLinks[csid][4]			# x of first connection point
			x1   = crossLinks[csid][5]			# x of second connection point

			xmt0 = xMT[mtl0][mt0][0]			# x of start point of first connection MT
			xac1 = xAc[acl1][ac1][0]			# x of start point of second connection Actin

			
			# Compute node ids of connection points
			nid0 = self.MTs[mtl0][mt0].n0 + np.round((x0-xmt0)/p.dl_MT).astype(int)
			nid1 = self.Acs[acl1][ac1].n0 + np.round((x1-xac1)/p.dl_Ac).astype(int)

			self.elements.append(CBARX2([self.nodes[nid1], self.nodes[nid0]],self.properties[4]));	
			self.elements[-1].localID = eCount;
			self.elements[-1].mechanism = self.mechanisms[0];
			self.elements[-1].mechanism.Initialize(self.elements[-1],self);
			eCount+=1;
		self.eidMTAcCL1 = eCount-1

		return (nCount, eCount)

