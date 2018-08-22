import os
from paraview.simple import *

# Compute filenames
main_folder = '/Users/rderooij/Documents/Stanford/GROWAN/Cython/SolidModels/BrainDevCirc_Q2/'
model_folder = 'RES_001/'

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
animationScene1.UpdateAnimationUsingDataTimeSteps()
animationScene1.PlayMode = 'Real Time'
animationScene1.Duration = 3

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.ResetCamera()

# add layout
viewLayout1 = GetLayout()
lineChartView1 = CreateView('XYChartView')
viewLayout1.AssignView(2, lineChartView1)


# Show density
densityDisplay = Show(paraview_00, renderView1)
ColorBy(densityDisplay, ('POINTS', 'Density'))
renderView1.ResetCamera()
densityDisplay.RescaleTransferFunctionToDataRange(True)
densityDisplay.SetScalarBarVisibility(renderView1, True)


# Plot density along x axis
plotOverLine1 = PlotOverLine(Input=paraview_00)
plotOverLine1.Source.Point1 = [0.0, 0.0, 0.0]
plotOverLine1.Source.Point2 = [60.00000000000000, 0.0, 0.0]
plotOverLine1Display = Show(plotOverLine1, lineChartView1)
plotOverLine1Display.CompositeDataSetIndex = [0]
plotOverLine1Display.UseIndexForXAxis = 0
plotOverLine1Display.XArrayName = 'arc_length'
plotOverLine1Display.SeriesLabel = ['Density', 'Density']
plotOverLine1Display.SeriesVisibility = ['Density']
lineChartView1.LeftAxisUseCustomRange = 1
lineChartView1.LeftAxisRangeMinimum = 0.0
lineChartView1.LeftAxisRangeMaximum = 0.5
lineChartView1.BottomAxisTitle = 'x [m]'
lineChartView1.BottomAxisUseCustomRange = 1
lineChartView1.BottomAxisRangeMinimum = 0.0
lineChartView1.BottomAxisRangeMaximum = 60.
lineChartView1.LeftAxisTitle = 'Density [kg/m3]'


"""
# create a new 'Threshold'
threshold1 = Threshold(Input=paraview_00)						# Create threshold
RenameSource('StateThreshold', threshold1)						# Rename
threshold1.Scalars = ['CELLS', 'State']							# Set parameter
threshold1.ThresholdRange = [0.0, 20.0]							# Set range
threshold1Display = Show(threshold1, renderView1)				# Show
ColorBy(threshold1Display, ('CELLS', 'State')) 					# Coloring
threshold1Display.RescaleTransferFunctionToDataRange(True) 		# Rescaling
threshold1Display.SetScalarBarVisibility(renderView1, False) 	# Hide scale bar
threshold1Display.LineWidth = 2.0								# Line thickness


# get color transfer function/color map for 'State'
stateLUT = GetColorTransferFunction('State')
stateLUT.RGBPoints = [1.0, 0.0, 0.0, 1.0, 20.0, 1.0, 0.0, 0.0]
stateLUT.ColorSpace = 'HSV'
stateLUT.NanColor = [0.498039, 0.498039, 0.498039]
stateLUT.ScalarRangeInitialized = 1.0
stateLUT.RescaleTransferFunction(0.0, 20.0)

# get opacity transfer function/opacity map for 'State'
statePWF = GetOpacityTransferFunction('State')
statePWF.Points = [1.0, 0.0, 0.5, 0.0, 20.0, 1.0, 0.5, 0.0]
statePWF.ScalarRangeInitialized = 1
statePWF.RescaleTransferFunction(0.0, 20.0)
"""


