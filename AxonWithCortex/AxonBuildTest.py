import numpy as np
from math import *			# for nan, pi
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os, sys
myPath = os.getcwd()
sys.path.insert(0, myPath+'/lsi')

from lsi import intersection
import params as p

def getGeometry():
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
				if(x1>p.lAxon):
					doContinue = False;
					x1 = np.minimum(x0+p.lMT0,p.lAxon)
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
	for i in range(1,int(p.lAxon/p.dlLink_MT_MT)):
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
				crossLinks.append([connCS[csID][0],id01,connCS[csID][1],id11,x0,x1])

				# Increment counter for cross links on these MT
				xMT[connCS[csID][0]][id01][2] +=1
				xMT[connCS[csID][1]][id11][2] +=1

		except:
			pass


	
	return {'xMT':xMT, 'yMT':yMT,'zMT':zMT,'crossLinks':crossLinks}

def getCortexGeometry():
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
			
			if x1>x0+1.e-6:
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

def getCL_MT_Ac(geomMT,geomAc):
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

	print connCS

	# Compute cross links
	crossLinks = [] # Every entry contains MTL1, MT1, AcL2, Ac2, x0, x1
	for i in range(1,int(p.lAxon/p.dlLink_MT_Ac)):
		x0 = i*p.dlLink_MT_Ac
		csID = np.random.randint(0,nCS)
		x1 = x0+connCS[csID][2]*np.tan(p.thLink_MT_Ac)

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


def plotGeometry(geomMT,geomAc,geomCL_Ac_Mt):
	xMT = geomMT['xMT']
	yMT = geomMT['yMT']
	zMT = geomMT['zMT']
	clMT = geomMT['crossLinks']

	xAc = geomAc['xAc']
	yAc  = geomAc['yAc']
	zAc = geomAc['zAc']
	clAc = geomAc['crossLinks']
	xAcRing = geomAc['xAcRing']
	nAcRing = len(xAcRing)

	clMTAc = geomCL_Ac_Mt['crossLinks']

	### =============================== ###
	###		MT and crossLinks in 3d 	###
	### ------------------------------- ###
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	colorMT = 'bgrcmy'
	colorAC = 'y'
	# MT
	for i in range(p.nMT):
		for j in range(len(xMT[i])):
			ax.plot([ xMT[i][j][0],xMT[i][j][1]],[yMT[i],yMT[i]],[zMT[i],zMT[i]],color=colorMT[j%len(colorMT)],linewidth=2.0)
	# CrossLinks
	for cl in clMT:
		ax.plot([ cl[4]     ,cl[5] ],\
				[ yMT[cl[0]],yMT[cl[2]] ],\
				[ zMT[cl[0]],zMT[cl[2]] ],'k',linewidth=0.5)
	
	### =============================== ###
	###		Actin and crossLinks in 3d 	###
	### ------------------------------- ###

	# Actin
	for i in range(p.nAcCirc):
		for j in range(len(xAc[i])):
			ax.plot([ xAc[i][j][0],xAc[i][j][1]],[yAc[i],yAc[i]],[zAc[i],zAc[i]],color=colorAC[j%len(colorAC)])
	# Actin rings
	thR = np.linspace(0,2*pi,100); yR = p.rAxon*np.cos(thR); zR = p.rAxon*np.sin(thR)
	for i in range(nAcRing):
		ax.plot(100*[xAcRing[i]],yR,zR,color=colorAC[i%len(colorAC)])

	# Crosslinks
	for cl in clAc:
		if cl[2]==-1:
			continue
		ax.plot([ cl[4]     ,cl[5] ],\
				[ yAc[cl[0]],yAc[cl[2]] ],\
				[ zAc[cl[0]],zAc[cl[2]] ],'k',linewidth=0.5)

	### =============================== ###
	###		Actin-MT  crossLinks in 3d 	###
	### ------------------------------- ###

	# Crosslinks
	for cl in clMTAc:
		ax.plot([ cl[4]     ,cl[5] ],\
				[ yMT[cl[0]],yAc[cl[2]] ],\
				[ zMT[cl[0]],zAc[cl[2]] ],'c',linewidth=0.5)

	"""
	for direction in (-1, 1):
		for point in np.diag(direction * p.lAxon * np.array([0,1,1])):
			ax.plot([point[0]], [point[1]], [point[2]], 'w')
	"""
	 
	ax.set_xlim(0.8*p.lAxon, 1.1 * p.lAxon)
	plt.savefig('axon1.eps',format='eps')
	ax.view_init(elev=0, azim=0)
	plt.savefig('axon2.eps',format='eps')
	ax.view_init(elev=0, azim=-90)
	plt.savefig('axon3.eps',format='eps')

	# Compute the number of cross-links for each MT element
	nCS_El = []
	nCS_0  = 0;		# Number of elements with 0 cross-links
	for i in xMT:
		for j in i:
			nCS_El.append(j[2])
			if(j[2]==0):
				nCS_0+=1
	print 'The axon model has '+ str(nCS_0) + ' unconnected elements'

	fig = plt.figure()
	plt.hist(nCS_El, 50, facecolor='green', alpha=0.75)
	plt.title((str(nCS_0)+' unconnected elements'))
	plt.xlabel('# of cross links at MT')
	plt.ylabel('frequency')
	plt.savefig('cs_hist.png',format='png')

# Main code
np.random.seed(0)

# geom contains:
# xMT 			list of lists. Outer lists contains all MT-lines, inner lists contain [x0,x1,nConn] of each MT element in that line
# yMT   		list of y coord of each MT-line
# zMT   		list of z coord of each MT-line
# crossLinks    list of [MT1_id, MT2_id, x0, x1] of each crossLink
geomMT = getGeometry()

# Compute cortex geometry, which contains:
# segments 	 list of segments with [[x0,y0,z0],[x1,y1,z1]] in each entry
# crossLinks list of cross links with [seg_id0, seg_id1, [x0,y0,z0],[x1,y1,z1]]
# xRing 	array of x-coordinates of all 'rings' of nodes (length = number of rings)
# idRing    array of ring id for every node (length = number of nodes)
geomAc = getCortexGeometry()

# Compute crosslinks between actin and MT, which contains:
geomCL_Ac_Mt = getCL_MT_Ac(geomMT,geomAc)



# Plot geometry
plotGeometry(geomMT,geomAc,geomCL_Ac_Mt)





