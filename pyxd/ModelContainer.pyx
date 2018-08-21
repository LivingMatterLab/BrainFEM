# -*- coding: utf-8 -*-
# ============================================ #
# This class contains the MODEL attributes 
# that do NOT change during the simulation,
# i.e. nodes, elements, materials, spc, etc.
# -------------------------------------------- #

cimport OutputHelper as oh

cdef class ModelContainer(object):
	def __init__(self):
		print "Initializing ModelContainer"
		self.materials = [] 
		self.properties = []  
		self.amplitudes = []  
		self.nodes = []  
		self.elements = [] 
		self.spc = [] 
		self.mpc = []
		self.loads = [] 
		self.contacts = [] 
		self.pressures = []
		self.mechanisms = []
		self.CreateOutputFolders()

	def __str__(self):
		oh.WriteToLog(self,"Model data:")
		oh.WriteToLog(self,"Number of materials:    " + str(self.nmat))
		oh.WriteToLog(self,"Number of properties:   " + str(self.nprop))
		oh.WriteToLog(self,"Number of amplitudes:   " + str(self.nampl))
		oh.WriteToLog(self,"Number of nodes:        " + str(self.nnode))
		oh.WriteToLog(self,"Number of elements:     " + str(self.nel))
		oh.WriteToLog(self,"Number of mpc:          " + str(self.nmpc))
		oh.WriteToLog(self,"Number of spc:          " + str(self.nspc))
		oh.WriteToLog(self,"Number of loads:        " + str(self.nloads))
		oh.WriteToLog(self,"Number of contacts:     " + str(self.ncontacts))
		oh.WriteToLog(self,"Number of pressures:    " + str(self.npressures))
		oh.WriteToLog(self,"Number of mechanisms:   " + str(self.nmech))
		oh.WriteToLog(self,"Number of dof:          " + str(self.ndof))
		oh.WriteToLog(self,"\n")

		return " "

	cpdef SetListLengths(self):
		self.nmat   = len(self.materials)
		self.nprop  = len(self.properties) 
		self.nampl  = len(self.amplitudes) 
		self.nnode  = len(self.nodes) 
		self.nel    = len(self.elements) 
		self.nspc   = len(self.spc)
		self.nmpc   = len(self.mpc)
		self.nloads = len(self.loads)
		self.ncontacts  = len(self.contacts)
		self.npressures = len(self.pressures)
		self.nmech  = len(self.mechanisms)
		self.ndof   = self.numDofPerNode*self.nnode

	cpdef BuildModel(self,object p):
		print "ERROR: BuildModel function not implemented in subclass"

	cpdef CreateOutputFolders(self):
		import os, glob, shutil
		cdef str folderMain
		cdef int REP_TOT

		#make save and work directories
		self.folderMain=os.getcwd()
		REP_TOT=1
		while os.path.exists(self.folderMain+'/RES_'+'{0:03d}'.format(REP_TOT))==True:
			REP_TOT=REP_TOT+1

		self.folderSimulation  = self.folderMain+'/RES_'+'{0:03d}'.format(REP_TOT)
		self.folderInput       = self.folderSimulation+'/Input'
		self.folderOutput      = self.folderSimulation+'/Output'
		self.folderPostProcess = self.folderSimulation+'/PostProcess'
		
		os.mkdir(self.folderSimulation)
		os.mkdir(self.folderInput)
		os.mkdir(self.folderOutput)
		os.mkdir(self.folderPostProcess)

		for f in glob.glob("*.py"):
			shutil.copy(f,self.folderInput+'/'+str(f))
		for f in glob.glob("*.pyx"):
			shutil.copy(f,self.folderInput+'/'+str(f))
		for f in glob.glob("*.pxd"):
			shutil.copy(f,self.folderInput+'/'+str(f))
		for f in glob.glob("*.xml"):
			shutil.copy(f,self.folderInput+'/'+str(f))

		# if ParaView.py exists, then write folder name in there
		try:
			f = open(self.folderInput + '/ParaView.py','r')
			filedata = f.read()
			f.close()
			f = open(self.folderInput + '/ParaView.py','w')
			f.write(filedata.replace("#REPLACE_HERE#","'RES_"+'{0:03d}'.format(REP_TOT)+"/'"))
			f.close()
		except:
			pass

	cpdef WriteStepOutput(self, Solver s):
		pass		# If not implemented for model, just pass