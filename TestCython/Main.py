# -*- coding: utf-8 -*-
import pyximport;
import os,sys; 

# Add paths
myPath = os.getcwd()
sys.path.insert(0, myPath+'/pyxd')

print myPath
pyximport.install()
from X import *


# Main code: create object of X class with some input
xObj = X(3)
