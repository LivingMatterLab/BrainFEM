import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import time

print "PostProcessing..."
startTime = time.time()
### ============================================== ###
### 					MT.txt 					   ### 
### ---------------------------------------------- ###

# Read "MT.txt" in the output folder
f = open('../Output/MT.txt', 'r')

stepTime  = []			# List of time points
stepLenMT = []			# List of total MT length for each time step
xMT = []				# List of list containing MT center x position
x1MT = []; x2MT=[];		# Lists of lists containing start and end positions for MT (to be removed later)

# Read file and write data into lists
line = f.readline()
countStep = 0;
while line:
	# New time step found:
	if line.startswith("Time"):
		# Read time point 
		stepTime.append(float(line.split()[1]))
		stepLenMT.append(0)

		# Skip line that says: "MT          X pos.      Y pos.      Z pos."
		f.readline()

		# Read x locations of all mt
		line = f.readline()

		countMT = 0;
		while not line=='\n':
			x1	 = float(line.split()[1]);
			x2   = float(line.split()[4]);
			stepLenMT[-1]+= (x2-x1)
			if countMT==0:
				xMT.append([0.5*(x1+x2)])
				x1MT.append([x1])
				x2MT.append([x2])
			else:
				xMT[countStep].append(0.5*(x1+x2))
				x1MT[countStep].append(x1)
				x2MT[countStep].append(x2)

			line = f.readline()

			countMT+=1;
		countStep+=1;
	line = f.readline()

# Convert xMT into array
xMT = np.array(xMT)
x1MT = np.array(x1MT)
x2MT = np.array(x2MT)

try:
	# Compute number of MT per cross section over entire length
	nVal = 1000;
	lVal = 1.1*np.max(xMT[-1,:])
	xMTpCS = np.linspace(0,lVal,nVal+1)
	dxVal = lVal/nVal
	nMTpCS = []

	for i in range(len(stepTime)):
		nMTpCS_i = np.zeros(nVal+1)
		for mt in range(len(x1MT[i])):
			x1 = x1MT[i][mt]
			x2 = x2MT[i][mt]
			nMTpCS_i[int(x1/dxVal):int(x2/dxVal)+1]+=1
		nMTpCS.append(nMTpCS_i)

	# Plot x vs time for each MT
	plt.figure()
	for i in range(len(xMT[0])):
		plt.plot(xMT[:,i],stepTime,'r',linewidth=1.0)
	#plt.axis('equal')
	#plt.axis('off')
	plt.xlabel('x')
	plt.ylabel('time')
	plt.gca().invert_yaxis()
	plt.savefig('../PostProcess/MT.eps',format='eps')
except:
	pass

### ============================================== ###
### 					FD.txt 					   ### 
### ---------------------------------------------- ###
stepTime = []; stepDisp = []; stepStretch = []; stepForce = []; stepPiola1 = []; stepNumMTCL = []; stepNumAcCL = []
f = open('../Output/FD.txt', 'r')

# Read file and write data into lists
line = f.readline()
line = f.readline();		# Skip first line
countStep = 0;
while line:
	lineSplit = line.split()
	stepTime.append(   float(lineSplit[1]))
	stepDisp.append(   float(lineSplit[2]))
	stepStretch.append(float(lineSplit[3]))
	stepForce.append(  float(lineSplit[4]))
	stepPiola1.append( float(lineSplit[5]))
	stepNumMTCL.append( float(lineSplit[6]))
	stepNumAcCL.append( float(lineSplit[7]))
	line = f.readline()


stepStretchRate = [];
stepViscosity   = [];
for i in range(len(stepTime)):

	i1 = np.maximum(0,i-int(0.05*len(stepTime)))
	i2 = np.minimum(i+int(0.05*len(stepTime)),len(stepTime)-1)

	try:
		m,b = np.polyfit(stepTime[i1:i2], stepStretch[i1:i2], 1)
	except:
		m = 1.;

	"""
	if (i==0 or i==len(stepTime)-1):
		stepStretchRate.append(np.nan);
	else:
		stepStretchRate.append((stepStretch[i+1]-stepStretch[i-1])/(stepTime[i+1]-stepTime[i-1]))
	"""
	stepStretchRate.append(m)
	stepViscosity.append(stepPiola1[i]/2./stepStretchRate[i]*stepStretch[i]**2)


fig, ax1 = plt.subplots()
ax1.plot(stepTime, stepStretch, 'b-')
ax1.set_xlabel('time')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('stretch', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')


ax2 = ax1.twinx()
ax2.plot(stepTime, stepForce, 'r-',markersize=10.0)
ax2.set_ylabel('force', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')

plt.savefig('../PostProcess/FD.eps',format='eps')

fig, ax1 = plt.subplots()
ax1.plot(stepTime, stepNumMTCL, 'b-')
ax1.set_xlabel('time')
ax1.set_ylabel('# of cross-links in MT', color='b')
ax1.set_ylim([0,1.2*ax1.get_ylim()[1]])
for tl in ax1.get_yticklabels():
    tl.set_color('b')

ax2 = ax1.twinx()
ax2.plot(stepTime, stepLenMT, 'r.',markersize=10.0)
ax2.set_ylabel('Total MT length', color='r')
ax2.set_ylim([0,1.2*ax2.get_ylim()[1]])
for tl in ax2.get_yticklabels():
    tl.set_color('r')

plt.savefig('../PostProcess/CL_MT.eps',format='eps')

fig, ax = plt.subplots()
ax.plot(stepTime, stepNumAcCL, 'b-')
ax.set_xlabel('time')
ax.set_ylabel('# of cross-links in Actin', color='b')
ax.set_ylim([0,1.2*ax.get_ylim()[1]])
for tl in ax.get_yticklabels():
    tl.set_color('b')
plt.savefig('../PostProcess/CL_Ac.eps',format='eps')

# Plot stretch rate and viscosity
fig, ax1 = plt.subplots()
ax1.plot(stepTime, stepStretchRate, 'r.',markersize=10.0)
ax1.set_xlabel('time')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('rate of stretch', color='r')
ax1.set_ylim([0,1.2*ax1.get_ylim()[1]])
for tl in ax1.get_yticklabels():
    tl.set_color('r')


ax2 = ax1.twinx()
ax2.plot([0,stepTime[-1]],[0.,0.],'k',linewidth=1.0)
ax2.plot(stepTime, stepViscosity, 'g.',markersize=10.0)
ax2.set_ylabel('viscosity', color='g')
for tl in ax2.get_yticklabels():
    tl.set_color('g')

plt.savefig('../PostProcess/Visc.eps',format='eps')



endTime = time.time()
print "Total time elapsed for post processing: "+str(endTime-startTime)
