# -*- coding: utf-8 -*-
import pyximport; 
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import time
import os, sys, traceback

# Add paths
myPath = os.getcwd()
sys.path.insert(0, myPath+'/../pyxd')
sys.path.insert(0, myPath+'/../pyxdX')
sys.path.insert(0, myPath+'/InputFiles')
pyximport.install(setup_args={"include_dirs":np.get_include()})

from DataContainer import *
from ModelContainer import *
from SOL2X import *
from DyneinModel_B3 import *
import ElementHelper as eh
import OutputHelper as oh
import params as p


#  ================================================  #
#  -----------  M A I N  P R O G R A M  -----------  #
#  ------------------------------------------------  #

for loadExt in [0.]:
	for irun in range(1):
		# Params for this run
		rseed = 0;#int(time.time()%1*10000000)			# Random seed
		p.loadExt = loadExt
		completed = 0.								# Simulation fraction completed

		try:
			# ModelContainer
			mc = DyneinModel_B3()
			np.random.seed(rseed)
			mc.BuildModel(p);
			print mc

			# DataContainer
			dc = DataContainer();
			dc.InitializeData(mc);

			# Plot
			#oh.ParaviewOutput(mc,dc) # Write paraview output

			# Solver
			sc = SOL2X(mc,dc,True,4)
			sc.maxIter = p.maxIter;
			sc.maxIterInc = p.maxIterInc;
			sc.maxStep = p.maxStep;
			sc.tEnd = p.tEnd;
			sc.dt0  = p.dt0;
			sc.dtMin = p.dtMin;
			sc.dtMax = p.dtMax;
			sc.tolerance = p.tolerance;
			sc.plotOutput = p.plotOutput;

			startTime = time.time()
			sc.Solve()
			endTime = time.time()
			oh.WriteToLog(mc,"Total time elapsed for solving: "+str(endTime-startTime))

			completed = dc.time/sc.tEnd;
		except:
			err = traceback.format_exc()		# Traceback of error
			oh.WriteToLog(mc,err)

	try:
		# Post process results
		os.chdir(mc.folderInput)
		execfile('PostProcess.py')
	except:
		err = traceback.format_exc()		# Traceback of error
		oh.WriteToLog(mc,err)

	
	strWrite = '{: <16s}'.format('Seed')       + '{: <16d}'.format(rseed) + '\n' +\
			   '{: <16s}'.format('loadExt')      + '{: <16.6f}'.format(p.loadExt) + '\n' +\
			   '{: <16s}'.format('Completed')  + '{: <16.6f}'.format(completed)
	oh.WriteToOutput(mc,'info.txt',strWrite)

	os.chdir(mc.folderMain)
