# -*- coding: utf-8 -*-
# ============================================ #
# This class contains DATA that changes during
# the simulation, i.e. dof, Fg, time, or data
# that is specific to type of solve (dt)
# Note that the latter should be implemented 
# in a solve class later.
# -------------------------------------------- #
cimport numpy as np
import numpy as np
import copy
from multiprocessing import Array

cdef class DataContainer(object):
	def __init__(self):
		self.Fg0      		= [];
		self.Fg       		= [];
		self.curvature0 	= [];
		self.curvature  	= [];
		self.theta0 		= [];
		self.theta  		= [];
		self.stretch0 		= [];
		self.stretch  		= [];
		self.nc0 			= [];
		self.nc  			= [];
		self.countFg0 		= [0];
		self.countFg  		= [0];
		self.countCurvature0= [0];
		self.countCurvature = [0];
		self.countTheta0    = [0];
		self.countTheta     = [0];
		self.countStretch0 	= [0];
		self.countStretch  	= [0];
		self.countNc0 		= [0];
		self.countNc  		= [0];
		self.OutputAnimation = [];

	cpdef InitializeData(self, ModelContainer mc):
		self.dof = Array('d', mc.ndof, lock=False)
		self.dof0 = Array('d', mc.ndof, lock=False)

		self.Rp = np.zeros(len(mc.spc)*mc.numDofPerNode)
		self.Rp0= np.zeros(len(mc.spc)*mc.numDofPerNode)

	cpdef ShareData(self):
		cdef int i,j,k,l,count

		###==========
		###    Fg
		###----------
		if self.countFg0[0]>0:
			self.Fg_  		= Array('d', self.countFg[0], lock=False)
			self.Fg0_ 		= Array('d', self.countFg0[0], lock=False)

			# Fill shared arrays with initializations
			count = 0;
			for i in range(len(self.Fg)):
				for j in range(len(self.Fg[i])):
					for k in range(len(self.Fg[i][0])):
						for l in range(len(self.Fg[i][1])):
							self.Fg_[count] = self.Fg[i][j][k,l];
							self.Fg0_[count]= self.Fg0[i][j][k,l];
							count = count+1;
		else:
			self.Fg_  		= np.zeros(0)
			self.Fg0_ 		= np.zeros(0)

		###==========
		###    Curvature and theta
		###----------
		if self.countCurvature0[0]>0:
			self.curvature_  = Array('d', self.countCurvature[0], lock=False)
			self.curvature0_ = Array('d', self.countCurvature0[0], lock=False)
			self.theta_  = Array('d', self.countTheta[0], lock=False)
			self.theta0_ = Array('d', self.countTheta0[0], lock=False)

			# Fill shared arrays with initializations
			count = 0;
			for i in range(len(self.curvature)):
				for j in range(len(self.curvature[i])):
					for k in range(len(self.curvature[i][0])):
						self.curvature_[count] = self.curvature[i][j][k];
						self.curvature0_[count]= self.curvature0[i][j][k];
						self.theta_[count] = self.theta[i][j][k];
						self.theta0_[count]= self.theta0[i][j][k];
						count = count+1;
		else:
			self.curvature_  = np.zeros(0)
			self.curvature0_ = np.zeros(0)
			self.theta_      = np.zeros(0)
			self.theta0_     = np.zeros(0)

		###==========
		###  stretch
		###----------
		if self.countStretch[0]>0:
			self.stretch_  	= Array('d', self.countStretch[0], lock=False)
			self.stretch0_ 	= Array('d', self.countStretch0[0], lock=False)

			count = 0;
			for i in range(len(self.stretch)):
				for j in range(len(self.stretch[i])):
					self.stretch_[count] = self.stretch[i][j];
					self.stretch0_[count]= self.stretch0[i][j];
					count = count+1;
		else:
			self.stretch_  		= np.zeros(0)
			self.stretch0_ 		= np.zeros(0)

		###==========
		###    nc
		###----------
		if self.countNc[0]>0:
			self.nc_  	= Array('d', self.countNc[0], lock=False)
			self.nc0_ 	= Array('d', self.countNc0[0], lock=False)

			count = 0;
			for i in range(len(self.nc)):
				for j in range(len(self.nc[i])):
					self.nc_[count] = self.nc[i][j];
					self.nc0_[count]= self.nc0[i][j];
					count = count+1;
		else:
			self.nc_  		= np.zeros(0)
			self.nc0_ 		= np.zeros(0)	
		
		# Delete initalizations because they are not used anymore
		del self.Fg0[:]
		del self.Fg[:]
		del self.curvature0[:]
		del self.curvature[:]
		del self.theta0[:]
		del self.theta[:]
		del self.stretch0[:]
		del self.stretch[:]
		del self.nc0[:]
		del self.nc[:]

	cpdef AddToData(self,listToAddTo, countToAddTo, dataToAdd):
		cdef int numDat, sizeDat, i

		# Add to list
		listToAddTo.append(dataToAdd);

		# Compute index of added data after listToAddTo will be transformed into 1d-array
		numDat = len(dataToAdd);
		try:
			sizeDat = dataToAdd[0].size
		except:
			# The data is just a single floating number
			sizeDat = 1;
		idData = [countToAddTo[0]+i*sizeDat for i in range(numDat)]

		# Update counter
		countToAddTo[0] += numDat*sizeDat

		return idData