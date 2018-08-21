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
from BarTest_B3 import *
import ElementHelper as eh
import OutputHelper as oh


#  ================================================  #
#  -----------  M A I N  P R O G R A M  -----------  #
#  ------------------------------------------------  #

# ModelContainer
mc = BarTest_B3()
mc.BuildModel(None);
print mc

# DataContainer
dc = DataContainer();
dc.InitializeData(mc);

#mc.TestElementMatrix(dc);


# Solver
sc = SOL2(mc,dc,False)
sc.maxIter = 10;
sc.maxIterInc = 6;
sc.maxStep = 2000;
sc.tEnd = 15.;
sc.dt0  = 1.;
sc.dtMin = 0.01;
sc.dtMax = 1.;


startTime = time.time()
sc.Solve()
endTime = time.time()
oh.WriteToLog(mc,"Total time elapsed for solving: "+str(endTime-startTime))

mc.PostProcess(dc)

print np.array(dc.stretch0_[:])
print np.array(dc.stretch_[:])
