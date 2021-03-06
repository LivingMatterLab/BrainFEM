# -*- coding: utf-8 -*-
import pyximport; 
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import time
import os, sys

# Add paths
myPath = os.getcwd()
sys.path.insert(0, myPath+'/../pyxd')
sys.path.insert(0, myPath+'/InputFiles')
pyximport.install(setup_args={"include_dirs":np.get_include()})

from DataContainer import *
from ModelContainer import *
from SOL2 import *
from BrainDevCirc_Q2 import *
import ElementHelper as eh
import OutputHelper as oh

import params as p


#  ================================================  #
#  -----------  M A I N  P R O G R A M  -----------  #
#  ------------------------------------------------  #
for E_inner in [1000.]:
	for E_ratio in [1.]:
		for kth_i in [47.]:
			for k_ratio in [1.]:
				for adv_alpha in [20.]:
					p.E_inner = E_inner
					p.E_ratio = E_ratio
					p.kth1_i  = kth_i
					p.kth2_i  = kth_i
					p.kth1_o  = kth_i*k_ratio
					p.kth2_o  = kth_i/k_ratio
					p.adv_R_alpha = adv_alpha

					try:
						# ModelContainer
						mc = BrainDevCirc_Q2()
						mc.BuildModel(p);
						print mc

						# DataContainer
						dc = DataContainer();
						dc.InitializeData(mc);

						# Solver
						sc = SOL2(mc,dc,True,4)
						sc.maxIter = p.maxIter;
						sc.maxIterInc = p.maxIterInc;
						sc.maxStep = p.maxStep;
						sc.tEnd = p.tEnd;
						sc.dt0  = p.dt0;
						sc.dtMin = p.dtMin;
						sc.dtMax = p.dtMax;
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
						# Postprocess
						from PostProcess import PostProcess
						PostProcess(mc)
					except:
						err = traceback.format_exc()		# Traceback of error
						oh.WriteToLog(mc,err)

					strWrite = '{: <16s}'.format('E_inner')     + '{: <16f}'.format(p.E_inner) + '\n' +\
							   '{: <16s}'.format('E_ratio')     + '{: <16.6f}'.format(p.E_ratio) + '\n' +\
							   '{: <16s}'.format('kth1_i')      + '{: <16.6f}'.format(p.kth1_i) + '\n' +\
							   '{: <16s}'.format('kth2_i')      + '{: <16.6f}'.format(p.kth2_i) + '\n' +\
							   '{: <16s}'.format('kth1_o')      + '{: <16.6f}'.format(p.kth1_o) + '\n' +\
							   '{: <16s}'.format('kth2_o')      + '{: <16.6f}'.format(p.kth2_o) + '\n' +\
							   '{: <16s}'.format('adv_R_alpha') + '{: <16.6f}'.format(p.adv_R_alpha) + '\n' +\
							   '{: <16s}'.format('Completed')  + '{: <16.6f}'.format(completed)
					oh.WriteToOutput(mc,'info.txt',strWrite)

					os.chdir(mc.folderMain)





