from ModelContainer cimport *
from Amplitude cimport *
from MAT_LE cimport * 
from MAT_VISC cimport * 
from PBAR cimport *
from PBAR1 cimport *
from NodeX cimport *
from CBARX cimport *
from SPC cimport *
from MPC cimport *
from LOAD cimport *
from MT cimport *
from MECH_EL01 cimport *
from MECH_EL02 cimport *
from MECH_MT02 cimport *
from Solver cimport *

cimport OutputHelper as oh
cimport ElementHelper as eh

from math import *
import numpy as np
cimport numpy as np

cdef class MTwithCortex_B3(ModelContainer):
	def __init__(self):
		super().__init__()
		print "Initializing MTwithCortex_B3"

	cpdef BuildModel(self,object p):
		self.numDofPerNode = 3
		self.lAxon = p.lAxon
		
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

		#  ================================================  #
		#    M A T E R I A L S  A N D  P R O P E R T I E S   #
		#  ------------------------------------------------  #
		# Create amplitude
		self.amplitudes.append(Amplitude(np.array([0.,p.tLoad,p.tEnd]),np.array([0.,1.,1.])));

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

		matProp['E'] = p.E_Med
		self.materials.append(MAT_LE(matProp))
		self.materials[-1].localID  = 14;				# Medium - elastic

		matProp = {'eta':p.eta_Ac}					
		self.materials.append(MAT_VISC(matProp))
		self.materials[-1].localID  = 15;				# Actin - viscous

		matProp = {'eta':p.eta_Med}
		self.materials.append(MAT_VISC(matProp))
		self.materials[-1].localID  = 16;				# Medium - viscous

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

		propProp = {'area':p.area_MT}
		self.properties.append(PBAR(self.materials[4],propProp))
		self.properties[-1].localID = 304; 				# Medium - elastic

		propProp['area'] = p.area_Ac;
		self.properties.append(PBAR(self.materials[5],propProp))
		self.properties[-1].localID = 205;				# Actin - viscous

		propProp = {'area':p.area_MT}
		self.properties.append(PBAR(self.materials[6],propProp))
		self.properties[-1].localID = 306;				# Medium - viscous

		for m in self.properties:
			print m

		# Create mechanism
		self.mechanisms.append(MECH_EL02(p.tCrea, p.tDest, p.maxInitStretch_Dyn, p.maxStretch_Dyn))
		self.mechanisms.append(MECH_MT02(p.tMTpoly,p.tMTstat,p.tMTdepoly,p.MTpolyRate,p.MTdepolyRate))
		
		#  ================================================  #
		#  ------  N O D E S  A N D  E L E M E N T S  -----  #
		#  ------------------------------------------------  #
		
		# Create nodes and MT elements
		self.MTs = [];				# List of list first and last nodeId for each MT, corresponding to xMT
		nCount = 0;						# Node count
		eCount = 0;						# Element count
		nCount, eCount, nMTViscNodes = self.CreateMT(p,geomMT, nCount, eCount)

		# Create nodes and elements for longitudinal actin filaments
		self.Acs = [];					# List of list first and last nodeId for each Actin, corresponding to xAc
		nCount, eCount, nAcViscNode, nidAcRing0, nidAcRing1 = self.CreateActin(p,geomAc, nCount, eCount)

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
		
		# Create dynein crossLink elements among MT
		nCount, eCount = self.CreateMTCrosslinks(p,geomMT,nCount,eCount);

		# Create myosin crosslink elements within the cortex
		nCount, eCount = self.CreateActinCrosslinks(p,geomAc,nidAcRing0, nidAcRing1, nCount,eCount);


		#  ================================================  #
		#  ----  B O U N D A R Y  C O N D I T I O N S  ----  #
		#  ------------------------------------------------  #
		# Compute nodeIDs of nodes at left (x=0) and right (x=lAxon) ends
		nidLeft = [];
		nidRight = [];
		for n in self.nodes:
			if np.abs(n.x)<1.e-6 and np.abs(n.y)<2*p.rAxon and np.abs(n.z)<2*p.rAxon:
				nidLeft.append(n.localID)
			if np.abs(n.x-p.lAxon)<1.e-6 and np.abs(n.y)<2*p.rAxon and np.abs(n.z)<2*p.rAxon:
				nidRight.append(n.localID)
		# Extract node that is loaded
		self.loadNode = self.nodes[nidRight[0]]

		# MPC
		# Constrain right nodes to move same amount in x direction
		for i in range(1,len(nidRight)):
			self.mpc.append(MPC([self.nodes[nidRight[i]],self.loadNode],[0,0],[-1.,1.]))

		# Constrain all actin rings to move same amount in x direction
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
		for nid in nMTViscNodes:
			self.spc.append(SPC(self.nodes[nid],[0],0.,self.amplitudes[0]))
		self.spc.append(SPC(self.nodes[nAcViscNode],[0],0.,self.amplitudes[0]))

		# Clamp all nodes in y,z direction
		for n in self.nodes:
			self.spc.append(SPC(n,[1,2],0.,self.amplitudes[0]))

		# Load master node at right in x direction
		self.loads.append(LOAD(self.loadNode,np.array([p.loadExt,0.,0.]),self.amplitudes[0]))
		
		# length of all matrices
		self.SetListLengths()
				
	cpdef getGeometry(self,object p):
		cdef np.ndarray[double,ndim=1] thMT, rMT, yMT, zMT
		cdef np.ndarray[np.int_t,ndim=1] lMT0
		cdef list xMT, connCS, crossLinks

		# Compute y,z of MicroTubules
		if p.nc==1:
			thMT = np.concatenate([np.array([0]),np.linspace(0.,2.*pi,7)])
			thMT = np.delete(thMT,[7])
			rMT  = np.concatenate([np.array([0]),p.rInner*np.ones(6)])

			print thMT
			print rMT
		elif p.nc==2:
			thMT = np.concatenate([np.array([0]),np.linspace(0.,2*pi,7),np.linspace(0.,2*pi,13)+0.*pi/12])
			thMT = np.delete(thMT,[7,20])
			rMT  = np.concatenate([np.array([0]),p.rInner*np.ones(6),np.array(6*[2.*p.rInner, sqrt(3.)*p.rInner])])
		else:
			raise "p.nc="+str(p.nc)+" is not supported!"

		yMT = rMT*np.cos(thMT)
		zMT = rMT*np.sin(thMT)

		# Compute possible cross-links
		connCS = []
		for i in range(p.nMT):
			for j in range(i+1,p.nMT):
				dist = np.sqrt((yMT[j]-yMT[i])**2+(zMT[j]-zMT[i])**2)
				if(dist<p.rConn):
					connCS.append([i,j,dist])
					connCS.append([j,i,dist])		# Same crosslink, but reversed order of nodes
		nCS = len(connCS)							# Number of potential crosslinks

		# Compute length of MT at x=0
		#lMT0 = (np.random.rand(p.nMT)*(p.lMT+p.lGap)).astype(int)
		#print lMT0
		lMT0 = np.linspace(0,p.lMTMax+p.lGap,p.nMT+1).astype(int)
		lMT0 = np.delete(lMT0,0)
		np.random.shuffle(lMT0)

		# For each MT-line, compute x0,x1 of each MT
		xMT = []
		for i in range(p.nMT):
			if(lMT0[i]<(p.lMTMax+p.lGap-p.lMT0)):
				x0 = lMT0[i]+p.lGap
				x1 = x0+p.lMTMax
			else:
				x0 = 0;
				x1 = lMT0[i]
			xMT.append([[x0,x1,0]]) 	# Note, the third entry is a counter for the number of crosslinks at this MT

			doContinue = True;
			while(doContinue):
				x0 = x1+p.lGap
				x1 = x0+p.lMTMax

				if(x0>p.lAxon):
					doContinue = False;
				else:
					if(x1>p.lAxon+p.dlGC):
						doContinue = False;
						x1 = np.minimum(x0+p.lMT0,p.lAxon+p.dlGC)
					xMT[i].append([x0,x1,0]) # Note, the third entry is a counter for the number of crosslinks at this MT

			# Correct length of first MT in this line if connected to the wall
			# This is because this first MT will never polymerize, so doesn't need additional
			# elements for polymerization
			if(xMT[i][0][0]==0):
				x1 = xMT[i][0][1]
				x1 = np.minimum(x1-p.lMTMax+p.lMT0,x1)
				xMT[i][0][1] = x1

		# Compute cross links
		crossLinks = [] # Every entry contains MTL1, MT1, MTL2, MT2, x0, x1
		for i in range(1,int((p.lAxon+p.dlGC)/p.dlLink_MT_MT)):
			x0 = i*p.dlLink_MT_MT
			csID = np.random.randint(0,nCS)
			x1 = x0+connCS[csID][2]*np.tan(p.thLink_MT_MT)

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
					# Check the fraction of both connections along its MT
					frac0 = (x0-xMT0[id01,0])/(xMT0[id01,1]-xMT0[id01,0])
					frac1 = (x1-xMT1[id11,0])/(xMT1[id11,1]-xMT1[id11,0])

					# For dynein motors, the walking side is towards the minus end of the MT
					# The walking side is the second node of the crosslink (node 1). 
					# This determines the order of nodes that define this crosslink
					if p.optionDynein==0 and frac0<frac1:
						# Yields extension
						crossLinks.append([connCS[csID][1],id11,connCS[csID][0],id01,x1,x0])
						"""
						# Order to yield extension
						if frac0>frac1:
							crossLinks.append([connCS[csID][0],id01,connCS[csID][1],id11,x0,x1])
						else:
							crossLinks.append([connCS[csID][1],id11,connCS[csID][0],id01,x1,x0])
						"""
					elif p.optionDynein==1:
						# No order
						crossLinks.append([connCS[csID][0],id01,connCS[csID][1],id11,x0,x1])
					elif p.optionDynein==2 and frac0>frac1:
						# Yields contraction
						crossLinks.append([connCS[csID][1],id11,connCS[csID][0],id01,x1,x0])

						"""
						# Order to yield contraction
						if frac0>frac1:
							crossLinks.append([connCS[csID][1],id11,connCS[csID][0],id01,x1,x0])
						else:
							crossLinks.append([connCS[csID][0],id01,connCS[csID][1],id11,x0,x1])
						"""
					else:
						# No element is created, so decrement counter
						# to compensate for increment on next line
						xMT[connCS[csID][0]][id01][2] -=1
						xMT[connCS[csID][1]][id11][2] -=1

					# Increment counter for cross links on these MT
					xMT[connCS[csID][0]][id01][2] +=1
					xMT[connCS[csID][1]][id11][2] +=1

			except:
				pass

		return {'xMT':xMT, 'yMT':yMT,'zMT':zMT,'crossLinks':crossLinks}

	cpdef getSimpleCortexGeometry(self, object p):
		cdef np.ndarray[double,ndim=1] yAc, zAc
		cdef list xAc, connCS, crossLinks

		# Compute y,z location of longitudinal actin segments
		yAc = np.array([p.rAxon*np.cos(0.)])
		zAc = np.array([p.rAxon*np.sin(0.)])
		xAc = [[[0,p.lAxon,0]]]
		connCS = []
		crossLinks = []
		xAcRing = [];

		return {'xAc':xAc, 'yAc':yAc,'zAc':zAc,'crossLinks':crossLinks,'xAcRing':xAcRing}

	cpdef getCortexGeometry(self, object p):
		cdef np.ndarray[double,ndim=1] thAc, yAc, zAc
		cdef list xAc, connCS, crossLinks

		# Compute y,z location of longitudinal actin segments
		thAc = np.linspace(0,2*pi,p.nAcCirc+1); thAc = np.delete(thAc,-1)
		yAc = p.rAxon*np.cos(thAc)
		zAc = p.rAxon*np.sin(thAc)
		# Every second Actin filament has larger radius
		for i in range(1,len(yAc),2):
			thAc[i] = thAc[i-1]
			yAc[i]  = yAc[i-1]*p.rFracActin
			zAc[i]  = zAc[i-1]*p.rFracActin

		# Compute x location of the actin rings
		xAcRing = [p.lActin/2.];  
		while xAcRing[-1]<p.lAxon:
			xAcRing.extend([xAcRing[-1]+p.dlActinRing])
		xAcRing = np.array(xAcRing)
		nAcRing = len(xAcRing)

		# For each Actin-line, compute x0,x1 of each Actin filament
		xAc = []
		for i in range(p.nAcCirc):
			for j in range(nAcRing/2):
				# x coordinates of actin ends
				x0 = np.maximum(0,       xAcRing[2*j+i%2]-p.lActin/2.)
				x1 = np.minimum(p.lAxon, xAcRing[2*j+i%2]+p.lActin/2.)
				
				if (x1>x0+1.e-6 and not (p.nAcCirc<=6 and (x1-x0)<p.lActin/4.)):
					xAc.append([[x0,x1,0]]) if j==0 else xAc[i].append([x0,x1,0])

		# Compute possible cross-links between longitudinal actin
		connCS = []
		dist = sqrt((yAc[0]-yAc[1])**2+(zAc[0]-zAc[1])**2)
		for i in range(0,p.nAcCirc,2):
			connCS.append([i,(i+1)%p.nAcCirc,dist])
			connCS.append([(i+1)%p.nAcCirc,i,dist])		# Same crosslink, but reversed order of nodes
		nCS = len(connCS)							# Number of potential crosslinks

		# Compute cross links between actin and rings
		crossLinks = [] # Every entry contains [AcL0, Ac0, AcL1, Ac1, x0, x1] OR [AcL0, Ac0, -1, Ring1, x0, th1]
		dx = p.lCL_Ac_Ac/sqrt(2.)
		dthInner = 2*np.arcsin(dx/p.rAxon)
		dthOuter = 2*np.arcsin(dx/p.rAxon/p.rFracActin)
		for i in range(p.nAcCirc):
			acl = xAc[i]
			for j in range(len(acl)):
				ac = acl[j]

				# Four cross links between actin filament and ring
				if xAcRing[2*j+i%2]<p.lAxon:
					if i%2 ==0:
						dth = dthInner
					else:
						dth = dthOuter
					crossLinks.append([i,j,-1,2*j+i%2,xAcRing[2*j+i%2]-dx,thAc[i]+dth])
					crossLinks.append([i,j,-1,2*j+i%2,xAcRing[2*j+i%2]+dx,thAc[i]+dth])
					crossLinks.append([i,j,-1,2*j+i%2,xAcRing[2*j+i%2]-dx,thAc[i]-dth])
					crossLinks.append([i,j,-1,2*j+i%2,xAcRing[2*j+i%2]+dx,thAc[i]-dth])


		# Compute crosslinks between longitudinal actin filaments
		for i in range(1,int(p.lAxon/p.dlLink_Ac_Ac)):
			x0 = i*p.dlLink_Ac_Ac
			csID = np.random.randint(0,nCS)
			x1 = x0+connCS[csID][2]*np.tan(p.thLink_Ac_Ac)

			# Check that there is no gap at either one of the MT
			xAc0 = np.asarray(xAc[connCS[csID][0]])
			xAc1 = np.asarray(xAc[connCS[csID][1]])

			try:
				id01 = [ n for n,i in enumerate(xAc0[:,0]) if i<=x0][-1]
				id02 = [ n for n,i in enumerate(xAc0[:,1]) if i>=x0][0]

				id11 = [ n for n,i in enumerate(xAc1[:,0]) if i<=x1][-1]
				id12 = [ n for n,i in enumerate(xAc1[:,1]) if i>=x1][0]

				# if id01==id02 and id11==id12, both endpoints of the crosslink are attached
				if(id01==id02 and id11==id12):

					
					# Check the fraction of both connections along its actin filament
					frac0 = (x0-xAc0[id01,0])/(xAc0[id01,1]-xAc0[id01,0])
					frac1 = (x1-xAc1[id11,0])/(xAc1[id11,1]-xAc1[id11,0])

					# Only create crosslink if both ends are pointing towards
					# center of MT, which is when myosin will create compression
					if( (frac0>0.5 and frac1<0.5 and x1>x0) or \
						(frac0<0.5 and frac1>0.5 and x1<x0)):
						crossLinks.append([connCS[csID][0],id01,connCS[csID][1],id11,x0,x1])

						# Increment counter for cross links on these MT
						xAc[connCS[csID][0]][id01][2] +=1
						xAc[connCS[csID][1]][id11][2] +=1
			except:
				pass


		if xAcRing[-1]>p.lAxon:
			xAcRing = np.delete(xAcRing,-1)
			nAcRing = len(xAcRing)

		return {'xAc':xAc, 'yAc':yAc,'zAc':zAc,'crossLinks':crossLinks,'xAcRing':xAcRing}

	cpdef getCL_MT_Ac(self, object p, dict geomMT, dict geomAc):
		cdef np.ndarray[double,ndim=1] yMT, zMT, yAc, zAc
		cdef list xMT, thMT, rMT, xAc, clMT, clAc, crossLinks, connCS
		cdef int i, j, nCS, id01, id02, id11, id12, csID
		cdef double dist, x0, x1

		# Extract variables
		xMT = geomMT['xMT']
		yMT = geomMT['yMT']
		zMT = geomMT['zMT']
		thMT = [np.arctan2(zMT[i],yMT[i])%(2*pi) for i in range(len(yMT))]
		rMT  = [sqrt(yMT[i]**2+zMT[i]**2) for i in range(len(yMT))]
		clMT = geomMT['crossLinks']

		xAc = geomAc['xAc']
		yAc  = geomAc['yAc']
		zAc = geomAc['zAc']
		clAc = geomAc['crossLinks']

		# Compute possible cross-links between MT and Actin
		connCS = []
		for i in range(p.nMT):
			for j in range(p.nAcCirc):
				dist = np.sqrt((yMT[i]-yAc[j])**2+(zMT[i]-zAc[j])**2)
				if(dist<p.rConn):
					connCS.append([i,j,dist])
		nCS = len(connCS)							# Number of potential crosslinks

		# Compute cross links
		crossLinks = [] # Every entry contains MTL1, MT1, AcL2, Ac2, x0, x1

		if nCS==0:
			return {'crossLinks':crossLinks}

		for i in range(1,int(p.lAxon/p.dlLink_MT_Ac)):
			x0 = i*p.dlLink_MT_Ac
			csID = np.random.randint(0,nCS)
			mpID = (1)**np.random.randint(0,2)	# Determines forward or backward tilt
			x1 = x0+mpID*connCS[csID][2]*np.tan(p.thLink_MT_Ac)

			# Check that there is no gap at either one of the MT
			xMT0 = np.asarray(xMT[connCS[csID][0]])
			xAc1 = np.asarray(xAc[connCS[csID][1]])

			try:
				id01 = [ n for n,i in enumerate(xMT0[:,0]) if i<=x0][-1]
				id02 = [ n for n,i in enumerate(xMT0[:,1]) if i>=x0][0]

				id11 = [ n for n,i in enumerate(xAc1[:,0]) if i<=x1][-1]
				id12 = [ n for n,i in enumerate(xAc1[:,1]) if i>=x1][0]

				# if id01==id02 and id11==id12, both endpoints of the crosslink are attached
				if(id01==id02 and id11==id12):
					crossLinks.append([connCS[csID][0],id01,connCS[csID][1],id11,x0,x1])
			except:
				pass

		return {'crossLinks':crossLinks}

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
				self.MTs[mtlID][-1].mechanism = self.mechanisms[1]
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
		
		# Create viscous elements connecting each MT to the clamped side
		initState = NoState
		initTimeToNextEvent = float("inf")
		for mtlID in range(len(yMT)):
			for mtID in range(len(xMT[mtlID])):
				# Create viscous bar element
				self.elements.append(CBARX([self.nodes[nViscNodes[mtlID]], self.nodes[self.MTs[mtlID][mtID].n0]], \
					                       self.properties[6],initState,initTimeToNextEvent));
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

			self.elements.append(CBARX([self.nodes[nid0], self.nodes[nid1]],self.properties[1]));	
			self.elements[-1].localID = eCount;
			self.elements[-1].mechanism = self.mechanisms[0];
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

				# Add segment as MT
				try:
					self.Acs[aclID].append(MT(nCount,nCount+nnodAc-1,nCount+nnodAc-1))
				except:
					self.Acs.append([MT(nCount,nCount+nnodAc-1,nCount+nnodAc-1)])
				if acID>0: 
					self.Acs[aclID][-1].mtMinus = self.Acs[aclID][-2]
					self.Acs[aclID][-2].mtPlus  = self.Acs[aclID][-1]
				self.Acs[aclID][-1].e0 = eCount
				self.Acs[aclID][-1].e2 = eCount+nelAc-1
				self.Acs[aclID][-1].e1 = self.Acs[aclID][-1].e2

				# Create elements
				for eid in range(nelAc):
					self.elements.append(CBARX([self.nodes[eid+nCount], self.nodes[eid+nCount+1]],self.properties[2]));
					self.elements[-1].state = Actin;
					self.elements[-1].localID = eid+eCount;
					self.nodes[eid+nCount].elPlus    = self.elements[-1]
					self.nodes[eid+nCount+1].elMinus = self.elements[-1]

					# If this is first element in this Actin filament, 
					# and this Actin filament is not the first along its line,
					# then assign corresponding elPlus and elMinus to nodes.
					# These elPlus and elMinus 'connect' two Actin filaments
					if(eid==0 and not acID==0):
						self.nodes[eid+nCount].elMinus = self.nodes[self.Acs[aclID][acID-1].n1].elMinus
						self.nodes[self.Acs[aclID][acID-1].n1].elPlus   = self.elements[-1]

				# Increment counts
				nCount += nnodAc;
				eCount += nelAc;

		# Create nodes and elements for circular actin filaments
		nR = np.ceil(2*pi*p.rAxon/p.dl_Ac).astype(int)
		thR = np.linspace(0,2*pi,nR+1); thR = np.delete(thR,-1)
		yR = p.rAxon*np.cos(thR); zR = p.rAxon*np.sin(thR)

		nidRing0 = []; nidRing1 = []
		for i in range(len(xAcRing)):
			# Create nodes
			for nid in range(nR):
				self.nodes.append(NodeX([xAcRing[i],yR[nid],zR[nid]]));
				self.nodes[-1].localID = nid+nCount;
				self.nodes[-1].dofID = range((nid+nCount)*self.numDofPerNode,(nid+nCount+1)*self.numDofPerNode);

			# Create elements
			for eid in range(nR):
				self.elements.append(CBARX([self.nodes[eid+nCount], self.nodes[(eid+1)%nR+nCount]],self.properties[2]));
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
			self.elements.append(CBARX([self.nodes[nid], self.nodes[nViscNode]],self.properties[6]));
			self.elements[-1].localID = eCount;
			eCount+=1

		# Attach end points of each longitudinal actin filaments to the viscNodes.
		for acl in self.Acs:
			for ac in acl:
				self.elements.append(CBARX([self.nodes[ac.n0], self.nodes[nViscNode]],self.properties[6]));
				self.elements[-1].localID = eCount;
				eCount+=1
	

		return (nCount, eCount, nViscNode, nidRing0, nidRing1)

	cpdef CreateActinCrosslinks(self,object p,dict geom, list nidRing0, list nidRing1, int nCount, int eCount):
		crossLinks = geom['crossLinks']

		nR = np.ceil(2*pi*p.rAxon/p.dl_Ac).astype(int)
		thR = np.linspace(0,2*pi,nR+1); thR = np.delete(thR,-1)
		yR = p.rAxon*np.cos(thR); zR = p.rAxon*np.sin(thR)
		
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
				

			self.elements.append(CBARX([self.nodes[nid0], self.nodes[nid1]],self.properties[3]));	
			self.elements[-1].localID = eCount;
			self.elements[-1].state = MotorAttachedTaut;
			eCount+=1;
		self.eidAcCL1 = eCount-1

		return (nCount, eCount)

	cpdef WriteStepOutput(self, Solver s):
		import params as p
		cdef str strWrite
		cdef double cLoad
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
		for i in range(self.eidMTCL0,self.eidMTCL1+1):
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

		# Count active actin crosslinks
		AcCLcount = 0
		for i in range(self.eidAcCL0,self.eidAcCL1+1):
			n0 = self.elements[i].nodes[0]

			# Check if this cross link is released
			if n0.localID == self.storageNodes[0].localID: continue

			AcCLcount +=1

		# Actin.txt contains location of all Actin filaments in this step
		strWrite = "Time " + '{: <15.6f}'.format(s.dc.time) +"\n"
		strWrite+= '{: <15}'.format('Actin') + '{: <15}'.format('N1.X') + '{: <15}'.format('N1.Y') + '{: <15}'.format('N1.Z')+ '{: <15}'.format('N2.X') + '{: <15}'.format('N2.Y') + '{: <15}'.format('N2.Z') +'\n'

		AcCount = 1
		for i in range(len(self.Acs)):
			for ac in self.Acs[i]:
				n0 = self.nodes[ac.n0]
				n1 = self.nodes[ac.n1]

				n0Dof = n0.Dof(s.dc)
				n1Dof = n1.Dof(s.dc)
				while len(n0Dof)<len(n0.loc): n0Dof.append(0.0);
				while len(n1Dof)<len(n1.loc): n1Dof.append(0.0);
				
				strWrite+= '{: <15d}'.format(AcCount) + \
				           '{: <15.6f}'.format(n0.x+n0Dof[0]) + \
				           '{: <15.6f}'.format(n0.y+n0Dof[1]) + \
				           '{: <15.6f}'.format(n0.z+n0Dof[2]) + \
				           '{: <15.6f}'.format(n1.x+n1Dof[0]) + \
				           '{: <15.6f}'.format(n1.y+n1Dof[1]) + \
				           '{: <15.6f}'.format(n1.z+n1Dof[2]) + \
				           '\n'

				AcCount +=1
		oh.WriteToOutput(self,'Actin.txt',strWrite)

		# FD.txt contains force-displacement at axon tip
		if s.dc.step==0:
			strWrite = '{: <16s}'.format('Step') + '{: <16s}'.format('Time') + \
					   '{: <16s}'.format('Displ. MT.') +  '{: <16s}'.format('Stretch MT.') + \
					   '{: <16s}'.format('Force') +  '{: <16s}'.format('Piola1')  + \
					   '{: <16s}'.format('# of MT CL') + \
					   '{: <16s}'.format('# of Ac CL') + '\n'
		else:
			strWrite = '';


		if p.optionBC==0:
			cLoad = self.amplitudes[0].Get(s.dc.time)*p.loadExt
		else:
			cLoad = -s.dc.Rp[-1]

		strWrite+= '{: <16d}'.format(s.dc.step) + \
		           '{: <16.6f}'.format(s.dc.time) + \
		           '{: <16.6e}'.format(self.growthConeNode.Dof(s.dc)[0]) + \
		           '{: <16.6f}'.format(1.+self.growthConeNode.Dof(s.dc)[0]/p.lAxon) + \
		           '{: <16.6e}'.format(cLoad) + \
		           '{: <16.6e}'.format(cLoad/pi/p.rAxon**2) + \
		           '{: <16d}'.format(CLcount) + \
		           '{: <16d}'.format(AcCLcount)

		oh.WriteToOutput(self,'FD.txt',strWrite)

	cpdef CalcLMT(self, Solver s):
		cdef double mtLength = 0.
		for mtl in self.MTs:
			for mt in mtl:
				n0 = self.nodes[mt.n0]
				n1 = self.nodes[mt.n1]
				mtLength+= n1.x-n0.x
				if not s==None:
					mtLength += n0.Dof(s.dc)[0]- n1.Dof(s.dc)[0]
		return mtLength

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
				mt.mechanism.Apply(mt,self,s)

		# Update total MT length
		self.currLMT = self.CalcLMT(s)

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
				eh.UpdateConnectivities(el,s)
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

			if mt.mtPlus:
				self.elements[mt.e1].nodes[1].elPlus = self.elements[mt.mtPlus.e0];
				self.nodes[mt.mtPlus.n0].elMinus = self.elements[mt.e1]
			else:
				self.elements[mt.e1].nodes[1].elPlus = None

		# Update total MT length
		self.currLMT = self.CalcLMT(s)

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
