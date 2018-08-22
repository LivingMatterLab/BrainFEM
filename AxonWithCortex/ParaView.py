import os
from paraview.simple import *

# Compute filenames
main_folder = '.../AxonWithCortex/'
model_folder = #REPLACE_HERE#

# Extract filename of .pvd file
pvdFile = None
folder = main_folder+model_folder+'PostProcess/'
for (dirpath, dirnames, filenames) in os.walk(folder):
	for file in filenames:
		if file.endswith('.pvd'):
			pvdFile = folder+file
	break

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
paraview_00 = PVDReader(FileName=pvdFile)
RenameSource(model_folder, paraview_00)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on animationScene1
animationScene1.PlayMode = 'Real Time'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ResetCamera()

# Threshold on State
threshold1 = Threshold(Input=paraview_00)						# Create threshold
RenameSource('StateThreshold', threshold1)						# Rename
threshold1.Scalars = ['CELLS', 'State']							# Set parameter
threshold1.ThresholdRange = [0.0, 20.0]							# Set range

# Threshold on Property - MT
threshold2 = Threshold(Input=threshold1)						# Create threshold
RenameSource('MT', threshold2)									# Rename
threshold2.Scalars = ['CELLS', 'PropID']						# Set parameter
threshold2.ThresholdRange = [50.0, 150.0]						# Set range
threshold2Display = Show(threshold2, renderView1)				# Show
ColorBy(threshold2Display, ('CELLS', 'State')) 					# Coloring
threshold2Display.RescaleTransferFunctionToDataRange(True) 		# Rescaling
threshold2Display.SetScalarBarVisibility(renderView1, False) 	# Hide scale bar
threshold2Display.LineWidth = 2.0								# Line thickness

# Threshold on Property - Cortex
threshold3 = Threshold(Input=threshold1)						# Create threshold
RenameSource('Cortex', threshold3)								# Rename
threshold3.Scalars = ['CELLS', 'PropID']						# Set parameter
threshold3.ThresholdRange = [150.0, 250.0]						# Set range
threshold3Display = Show(threshold3, renderView1)				# Show
ColorBy(threshold3Display, ('CELLS', 'State')) 					# Coloring
threshold3Display.RescaleTransferFunctionToDataRange(True) 		# Rescaling
threshold3Display.SetScalarBarVisibility(renderView1, False) 	# Hide scale bar
threshold3Display.LineWidth = 2.0								# Line thickness

# Threshold on Property - MT-Cortex interface
threshold4 = Threshold(Input=threshold1)						# Create threshold
RenameSource('MTCortInterface', threshold4)								# Rename
threshold4.Scalars = ['CELLS', 'PropID']						# Set parameter
threshold4.ThresholdRange = [250.0, 350.0]						# Set range
threshold4Display = Show(threshold4, renderView1)				# Show
ColorBy(threshold4Display, ('CELLS', 'State')) 					# Coloring
threshold4Display.RescaleTransferFunctionToDataRange(True) 		# Rescaling
threshold4Display.SetScalarBarVisibility(renderView1, False) 	# Hide scale bar
threshold4Display.LineWidth = 2.0								# Line thickness

renderView1.ResetCamera()

# get color transfer function/color map for 'State'
stateLUT = GetColorTransferFunction('State')
f = open(main_folder+model_folder+'Input/ParaViewColors.xml', "r"); cdata = f.read()
stateLUT.ApplyColorMap(cdata)
