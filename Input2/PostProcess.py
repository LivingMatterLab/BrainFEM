import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import time
from math import *

def PostProcess(model):
	print "PostProcessing..."
	startTime = time.time()

	densThreshold = 1.5*0.04
	grayColor =  [130./255.,  80./255.,  45./255.]
	whiteColor = [240./255., 230./255., 220./255.]

	# Experimental prediction of brain volume
	aExp =   2.65
	bExp = -15.15
	t0Exp = 10.

	tEnd = 19.

	### ============================================== ###
	### 				READ vol.txt	 			   ### 
	### ---------------------------------------------- ###
	f = open(model.folderOutput+'/vol.txt', 'r')
	stepTime = []; stepVolume = [];stepGrayVolume = []; stepWhiteVolume = [];  

	# Read file and write data into lists
	line = f.readline()
	countStep = 0;
	while line:
		# New time step found:
		if line.startswith("Time"):
			# Read time point 
			stepTime.append(float(line.split()[1]))
			stepVolume.append(0)
			stepGrayVolume.append(0)
			stepWhiteVolume.append(0)

			# Skip line that says: "Element        Volume         Density   "
			f.readline()

			# Read volume of all elements
			line = f.readline()

			while not line=='\n':
				elVol = float(line.split()[1])
				elDens = float(line.split()[2])

				stepVolume[-1]+=elVol;
				if elDens>densThreshold:
					stepGrayVolume[-1]+=elVol;
				else:
					stepWhiteVolume[-1]+=elVol;

				line = f.readline()
		line = f.readline()

	### ============================================== ###
	### 			Compute Volume Error 	 	   	   ### 
	### ---------------------------------------------- ###
	vol3dExp  = (aExp*(np.array(stepTime)+t0Exp)+bExp)**3.
	vol2dExp = (vol3dExp/2./pi)**(2./3)*pi/2.
	timeWeight = 0.5*(1.-np.array(stepTime)/tEnd)+0.5
	volumeDiff = np.multiply(timeWeight, vol2dExp-np.array(stepVolume))
	volumeError = np.linalg.norm(volumeDiff)
	# Penalize non-completed simulations
	volumeError/= (stepTime[-1]/tEnd)**5.

	# Penalize absence of grayVolume
	grayRatio = np.array(stepGrayVolume)/np.array(stepVolume)
	volumeError += 10000*max([0.,0.3-grayRatio[-1]])**2

	# Penalize late onset of gray development
	grayRatioLimit = 0.05
	idg = np.where(grayRatio<grayRatioLimit)[0]
	if len(idg)<=1 or idg[-1]==len(stepTime)-1:
		tGray = tEnd
	else:
		dRat = grayRatio[idg[-1]+1]-grayRatio[idg[-1]]
		dTim = stepTime[idg[-1]+1]-stepTime[idg[-1]]
		tGray = stepTime[idg[-1]]+dTim/dRat*(grayRatioLimit-grayRatio[idg[-1]])
	volumeError += 50*(tGray/tEnd)**2




	### ============================================== ###
	### 			Experimental prediction 	 	   ### 
	### ---------------------------------------------- ###
	tExp = np.linspace(0,max(stepTime),1000)
	volExp  = (aExp*(tExp+t0Exp)+bExp)**3.
	areaExp = (volExp/2./pi)**(2./3)*pi/2.

	# Plot volume over time
	fig, ax1 = plt.subplots()
	ax1.fill_between(stepTime, np.zeros(len(stepVolume)),stepWhiteVolume, color=whiteColor)
	ax1.fill_between(stepTime, stepWhiteVolume,stepVolume, color=grayColor)
	ax1.plot(stepTime, stepWhiteVolume, 'k')
	ax1.plot(stepTime, stepVolume, 'k')
	ax1.plot(tExp, areaExp, 'b',lw=3)
	ax1.set_xlabel('Time [s]')
	ax1.set_ylabel('Volume [m3]')
	ax1.set_ylim(0,ax1.get_ylim()[1])
	plt.savefig(model.folderPostProcess+'/vol.eps',format='eps')

	### ============================================== ###
	### 				Initialize Axes		   		   ### 
	### ---------------------------------------------- ###
	# Two subplots, unpack the axes array immediately
	fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,1, sharex=True)
	lineColor = [255./255,205/255.,153/255.]
	axisColor = [128/255., 0/255.,0/255.]

	fig.set_size_inches(6,6)

	# Axis 1: Stiffness E
	ax1.set_ylabel(r'$E$ []', multialignment='center',color=axisColor,labelpad=8)

	# Axis 2: Diffusivity D
	ax2.set_ylabel(r'$D$ []', multialignment='center',color=axisColor,labelpad=16)

	# Axis 3: Source Grho
	ax3.set_ylabel(r'$f^{\rho}$ []', multialignment='center',color=axisColor,labelpad=8)

	# Axis 4: Advection speed
	ax4.set_ylabel(r'$v_{adv}$ []', multialignment='center',color=axisColor,labelpad=16)

	# Axis 5: kth1, kth2
	ax5.set_ylabel(r'$k_{th}$ []', multialignment='center',color=axisColor,labelpad=16)
	ax5.set_xlabel(r'Radius',color=axisColor)

	for ax in [ax1,ax2,ax3,ax4,ax5]:
		ax.spines['left'].set_color(axisColor)
		ax.spines['right'].set_color(axisColor)
		ax.spines['bottom'].set_color(axisColor)
		ax.spines['top'].set_color(axisColor)
		for tl in ax.get_yticklabels():
		    tl.set_color(axisColor)
		for tl in ax.get_xticklabels():
		    tl.set_color(axisColor)


	### ============================================== ###
	### 			Plot for each element		   	   ### 
	### ---------------------------------------------- ###
	for el in model.elements:
		# Compute center of element
		locN = np.array([el.nodes[0].loc, el.nodes[1].loc,\
		                 el.nodes[2].loc, el.nodes[3].loc])

		locAvg = np.array([np.mean(locN[:,i]) for i in range(len(locN[0]))])
		xA = locAvg[0]
		yA = locAvg[1]
		rA = np.sqrt(xA**2+yA**2)
		if yA<rA/10.:
			continue


		ax1.plot(rA,el.property.material.E,'.',color = lineColor,markersize=10.0)
		ax2.plot(rA,el.property.material.D,'.',color = lineColor,markersize=10.0)
		ax3.plot(rA,el.property.Grho,'.',color = lineColor,markersize=10.0)
		ax4.plot(rA,el.property.advSpeed,'.',color = lineColor,markersize=10.0)
		ax5.plot(rA,el.property.kth1,'.',color = lineColor,markersize=10.0)
		ax5.plot(rA,el.property.kth2,'.',color = axisColor,markersize=10.0)



	### ============================================== ###
	### 				Save Figure		   		   	   ### 
	### ---------------------------------------------- ###
	plt.savefig(model.folderInput+'/InputParam.eps',format='eps')
	plt.close()

	endTime = time.time()
	print "Total time elapsed for post processing: "+str(endTime-startTime)

	return volumeError
