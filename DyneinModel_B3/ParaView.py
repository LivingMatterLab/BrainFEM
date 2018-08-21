import os
from paraview.simple import *

# Compute filenames
main_folder = '/Users/rderooij/Documents/Stanford/GROWAN/Cython/StudentCodes/IngridCode/MainCode/DyneinModel_B3/'
model_folder = #REPLACE_HERE#

paraFiles = []
folder = main_folder+model_folder+'PostProcess/'
for (dirpath, dirnames, filenames) in os.walk(folder):
	for file in filenames:
		if file.endswith('.vtu'):
			paraFiles.append(folder+file)
	break


#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
paraview_00 = XMLUnstructuredGridReader(FileName=paraFiles)
RenameSource(model_folder, paraview_00)
paraview_00.CellArrayStatus = ['Jg', 'PropID', 'State']
paraview_00.PointArrayStatus = ['Displacement']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ResetCamera()

# create a new 'Threshold'
threshold1 = Threshold(Input=paraview_00)						# Create threshold
RenameSource('StateThreshold', threshold1)						# Rename
threshold1.Scalars = ['CELLS', 'State']							# Set parameter
threshold1.ThresholdRange = [-1.0, 10.0]							# Set range
threshold1Display = Show(threshold1, renderView1)				# Show
ColorBy(threshold1Display, ('CELLS', 'State')) 					# Coloring
threshold1Display.RescaleTransferFunctionToDataRange(True) 		# Rescaling
threshold1Display.SetScalarBarVisibility(renderView1, False) 	# Hide scale bar
threshold1Display.LineWidth = 2.0								# Line thickness
renderView1.ResetCamera()

# get color transfer function/color map for 'State'
stateLUT = GetColorTransferFunction('State')
stateLUT.RGBPoints = [1.0, 0.0, 0.0, 1.0, 10.0, 1.0, 0.0, 0.0]
stateLUT.ColorSpace = 'HSV'
stateLUT.NanColor = [0.498039, 0.498039, 0.498039]
stateLUT.ScalarRangeInitialized = 1.0
stateLUT.RescaleTransferFunction(0.0, 10.0)

# get opacity transfer function/opacity map for 'State'
statePWF = GetOpacityTransferFunction('State')
statePWF.Points = [1.0, 0.0, 0.5, 0.0, 10.0, 1.0, 0.5, 0.0]
statePWF.ScalarRangeInitialized = 1
statePWF.RescaleTransferFunction(0.0, 10.0)



