from ModelContainer cimport *
from Amplitude cimport *
from MAT_NEOH cimport * 
from PSOLID71 cimport *
from Node cimport *
from CQUAD7 cimport *
from SPC cimport *
from LOAD cimport *
from Solver cimport *
import OutputHelper as oh

from math import *
import numpy as np
cimport numpy as np

cdef class BrainDevCirc_Q2(ModelContainer):
	def __init__(self):
		super().__init__()
		print "Initializing BrainDevCirc_Q2"

	cpdef BuildModel(self, object p):
		self.numDofPerNode = 3
		self.numDispDofPerNode = 2


		# Read INPUT
		cdef str filename = 'nc_'+str(p.nelc)+'.inp'
		nod, elem, cortex_elem, subcor_elem, core_elem = self.ReadAbaqusFile(filename)
		nel = len(elem)
		radius = np.max(nod[:,1])

		
		# Amplitudes
		r_amp = np.linspace(0.,radius,1000)
		t_amp = np.linspace(0.,p.tEnd,1000)
		rho_amp = np.linspace(0.,.1,1000);

		src_R  =  1-  np.exp(p.src_R_alpha*(r_amp-p.src_R_center*radius))/ \
		         (1+ (np.exp(p.src_R_alpha*(r_amp-p.src_R_center*radius))) );

		src_T  =  1-  np.exp(p.src_T_alpha*(t_amp-p.src_T_center))/ \
		         (1+ (np.exp(p.src_T_alpha*(t_amp-p.src_T_center))) );

		adv_R  =  1-  np.exp(p.adv_R_alpha*(r_amp-p.adv_R_center*radius))/ \
		         (1+ (np.exp(p.adv_R_alpha*(r_amp-p.adv_R_center*radius))) );

		adv_RHO=      np.exp(p.adv_RHO_alpha*(rho_amp-p.adv_RHO_center))/ \
		         (1+ (np.exp(p.adv_RHO_alpha*(rho_amp-p.adv_RHO_center))) );

		d_adv_RHO=  p.adv_RHO_alpha * \
		              np.exp(p.adv_RHO_alpha*(rho_amp-p.adv_RHO_center))/ \
		         (1+ (np.exp(p.adv_RHO_alpha*(rho_amp-p.adv_RHO_center))) )**2;         

		self.amplitudes.append(Amplitude(np.array([0.,p.tLoad,p.tEnd]),np.array([0.,1.,1.])));
		self.amplitudes.append(Amplitude(r_amp,src_R));
		self.amplitudes.append(Amplitude(t_amp,src_T));
		self.amplitudes.append(Amplitude(r_amp,adv_R));
		self.amplitudes.append(Amplitude(rho_amp,adv_RHO,d_adv_RHO));

		# Materials
		cdef dict matProp = {'E':p.E_inner,'nu':0.0,'D':0}
		for pid in range(p.nMatProp):
			ampl3 = self.amplitudes[3].Get(1.*(pid+1)/p.nMatProp*radius)

			matProp['E'] = p.E_inner*p.E_ratio*(1.-(1.-1./p.E_ratio)*ampl3)
			matProp['D'] = p.D*(1-0.95*ampl3)
			self.materials.append(MAT_NEOH(matProp))
			self.materials[-1].localID  = 10+pid;

		# Property
		cdef dict propProp = {'kth1':p.kth1_i,'kth2':p.kth2_i,'alpha':p.alpha,\
		                      'Grho':p.Grho,'advSpeed':p.advSpeed,\
		                      'amplRhoSource':self.amplitudes[2],\
		                      'amplAdvThreshold':self.amplitudes[4]}
		for pid in range(p.nMatProp):
			ampl1 = self.amplitudes[1].Get(1.*(pid+1)/p.nMatProp*radius)
			ampl3 = self.amplitudes[3].Get(1.*(pid+1)/p.nMatProp*radius)

			propProp['kth1'] = p.kth1_o*(1.-(1.-p.kth1_i/p.kth1_o)*ampl3)
			propProp['kth2'] = p.kth2_o*(1.-(1.-p.kth2_i/p.kth2_o)*ampl3)	                         
			propProp['Grho'] = p.Grho * ampl1
			propProp['advSpeed'] = p.advSpeed * ampl3
			self.properties.append(PSOLID71(self.materials[pid],propProp))
			self.properties[-1].localID = 100+pid;

		# Nodes
		cdef int nid;
		cdef double xi, yi
		for nid in range(len(nod)):
			xi = nod[nid,1]
			yi = nod[nid,2]
			thi = np.arctan2(yi,xi);

			if(np.abs(thi-pi/2.)<=p.pert_width and yi>=p.subcort_frac*radius-1.e-6):
				yi += p.pert*(1.-p.subcort_frac)*radius

			self.nodes.append(Node([xi,yi]));
			self.nodes[nid].localID = nid;
			self.nodes[nid].dofID = range(nid*self.numDofPerNode,(nid+1)*self.numDofPerNode);


		# Elements
		cdef int count, eid
		cdef np.ndarray[np.int_t, ndim=1] elnod
		cdef np.ndarray[double, ndim=1] locAvg, elemDir
		cdef np.ndarray[double, ndim=2] locN
		cdef double xA, yA
		count = 0

		for eid in range(nel):
			elnod = elem[eid,1:5]

			# Compute center of element
			locN = np.array([self.nodes[elnod[0]].loc, self.nodes[elnod[1]].loc,\
			                 self.nodes[elnod[2]].loc, self.nodes[elnod[3]].loc])

			locAvg = np.array([np.mean(locN[:,i]) for i in range(len(locN[0]))])
			xA = locAvg[0]
			yA = locAvg[1]

			# Element direction
			elemDir = self.GetElemDir(xA,yA,p.elemDir)

			# Property id
			pid = int(np.round(sqrt(xA**2+yA**2)/radius*p.nMatProp))

			# Create element
			self.elements.append(CQUAD7([self.nodes[int(ii)] for ii in elnod],elemDir,self.properties[pid]));
			self.elements[count].localID = count;

			count+=1

		
		# SPC
		# Clamp left nodes in y
		nidC = -1;
		xMin = 1e9;
		for n in self.nodes:
			if abs(n.y)<1.e-6:
				self.spc.append(SPC(n,[1],0.,self.amplitudes[0]))
				if abs(n.x)<xMin:
					nidC = n.localID
					xMin = abs(n.x)
		# Clamp bottom center node in x
		self.spc.append(SPC(self.nodes[nidC],[0],0.,self.amplitudes[0]))
		print self.nodes[nidC]

		# LOADS
		# - 
		
		# length of all matrices
		self.SetListLengths()

	cpdef GetElemDir(self, double x, double y, str option):
		cdef np.ndarray[double,ndim=1] elemDir
		
		if(option=='random'):
			th = 2*pi*np.random.rand()
			elemDir = np.array([cos(th), sin(th)])
		elif(option=='horizontal'):
			elemDir = np.array([1.,0.])
		elif(option=='vertical'):
			elemDir = np.array([0.,1.])
		elif(option=='diagonal'):
			elemDir = np.array([1.,1.])/sqrt(2.)
		elif(option=='radial'):
			elemDir = np.array([x,y])
			elemDir/= np.linalg.norm(elemDir)
		else:
			print "ERROR:  elemDir option ", option, " is not supported"
		
		return elemDir

	cpdef ReadAbaqusFile(self,str filename):
		# Open file
		fd = open('InputFiles/'+filename,'r')

		# Maximum number of lines
		maxLine = 1000000;

		##########
		### NODES
		##########
		fd.seek(0)
		iterr = 0
		found = 0
		line = fd.readline()
		while line!="" and iterr<maxLine:
			iterr+=1
			if '*Node' in line:
				found = 1
				break
			line = fd.readline()

		#Read lines
		iterr = 0
		line = fd.readline()
		nodes = []
		while line!="" and iterr<maxLine:
			data = line.strip().split(',')
			try:
				nodes.append([float(data[i]) for i in range(len(data))]);
			except:
				break


			line = fd.readline()
			iterr+=1

		# Write as array
		nodes=np.array([np.array(nodi) for nodi in nodes])	

		##########
		# ELEMENTS
		##########
		fd.seek(0)
		iterr = 0
		found = 0
		line = fd.readline()
		while line!="" and iterr<maxLine:
			iterr+=1
			if '*Element' in line:
				found = 1
				break
			line = fd.readline()

		#Read lines
		iterr = 0
		line = fd.readline()
		elements = []
		while line!="" and iterr<maxLine:
			data = line.strip().split(',')
			try:
				elements.append([int(data[i]) for i in range(len(data))]);
			except:
				break

			line = fd.readline()
			iterr+=1

		# Write as array
		elements=np.array([np.array(eli) for eli in elements])

		##########
		# CORTEX ELEMENTS
		##########
		fd.seek(0)
		iterr = 0
		line = fd.readline()
		while line!="" and iterr<maxLine:
			iterr+=1
			if '*Elset, elset=ELLIPSE_FACE-CORTEX, generate' in line:
				line = fd.readline()
				try:
					data = line.strip().split(',')
					cortex_elements = range(int(data[0]),int(data[1])+1,int(data[2]))
				except:
					break
				break
			line = fd.readline()
			iterr+=1

		# Write as array
		cortex_elements=np.array(cortex_elements)

		##########
		# SUBCORTEX ELEMENTS
		##########
		fd.seek(0)
		iterr = 0
		line = fd.readline()
		while line!="" and iterr<maxLine:
			iterr+=1
			if '*Elset, elset=ELLIPSE_FACE-SUBCORTEX, generate' in line:
				line = fd.readline()
				try:
					data = line.strip().split(',')
					subcortex_elements = range(int(data[0]),int(data[1])+1,int(data[2]))
				except:
					break
				break
			line = fd.readline()
			iterr+=1

		# Write as array
		subcortex_elements=np.array(subcortex_elements)

		##########
		# CORE ELEMENTS
		##########
		fd.seek(0)
		iterr = 0
		line = fd.readline()
		while line!="" and iterr<maxLine:
			iterr+=1
			if '*Elset, elset=ELLIPSE_FACE-CORE, generate' in line:
				line = fd.readline()
				try:
					data = line.strip().split(',')
					core_elements = range(int(data[0]),int(data[1])+1,int(data[2]))
				except:
					break
				break
			line = fd.readline()
			iterr+=1

		# Write as array
		core_elements=np.array(core_elements)

		##########
		# FINAL CONSIDERATIONS
		##########
		fd.close()

		# Account for the fact that abaqus starts at 1, python at 0
		for i in range(len(nodes)):
			nodes[i,0]-=1 
		
		elements-=1
		cortex_elements-=1
		subcortex_elements-=1
		core_elements-=1

		return (nodes, elements, cortex_elements, subcortex_elements,core_elements)

	cpdef WriteStepOutput(self, Solver s):
		cdef str strWrite
		cdef Element el;
		cdef int i, j

		# vol.txt contains volume and density in all elements in this step
		strWrite = "Time " + '{: <15.6f}'.format(s.dc.time) +"\n"
		strWrite+= '{: <15}'.format('Element') + '{: <15}'.format('Volume') + '{: <15}'.format('Density') + '\n'

		for el in self.elements:
			strWrite+= '{: <15d}'.format(el.localID) + \
			           '{: <15.6f}'.format(el.getVolume(s.dc)) + \
			           '{: <15.6f}'.format(el.getDensity(s.dc)) + \
			           '\n'
		oh.WriteToOutput(self,'vol.txt',strWrite)
