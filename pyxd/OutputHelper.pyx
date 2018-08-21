# -*- coding: utf-8 -*-
# ============================================ #
# This file contains general FUNCTIONS used to
# plot meshed or create output files
# -------------------------------------------- #
cimport numpy as np
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from ModelContainer cimport *
from DataContainer cimport *
from Element cimport *
cimport ElementHelper as eh
import vtk
import os


#Generate Paraview file
cpdef ParaviewOutput(ModelContainer mc, DataContainer dc, int it=0):
	cdef list nodDof

	# Initialize output variables
	Jg = vtk.vtkDoubleArray()
	Jg.SetName("Jg")
	Jg.SetNumberOfTuples(mc.nel)

	Stretch = vtk.vtkDoubleArray()
	Stretch.SetName("Stretch")
	Stretch.SetNumberOfTuples(mc.nel)

	PropID = vtk.vtkIntArray()
	PropID.SetName("PropID")
	PropID.SetNumberOfTuples(mc.nel)

	State = vtk.vtkIntArray()
	State.SetName("State")
	State.SetNumberOfTuples(mc.nel)
	
	uvw = vtk.vtkDoubleArray()
	uvw.SetName("Displacement")
	uvw.SetNumberOfTuples(mc.nnode*3)
	uvw.SetNumberOfComponents(3)

	rho = vtk.vtkDoubleArray()
	rho.SetName("Density")
	rho.SetNumberOfTuples(mc.nnode)

	xyz = vtk.vtkPoints()
	xyz.SetDataTypeToDouble()
	xyz.Allocate(mc.nnode)

	# Extract current locations and displacements
	for nod in mc.nodes:
		nodDof = nod.Dof(dc)
		# Fill nodDof with zeros if necessary (only when len(nodDof)<len(nod.loc))
		while (len(nodDof)<len(nod.loc)):
			nodDof.append(0.0);

		if len(nod.loc)==2:	
			xyz.InsertPoint(nod.localID, nod.loc[0]+nodDof[0],  nod.loc[1]+nodDof[1], 0.)
			uvw.SetTuple3(nod.localID, nodDof[0], nodDof[1], 0.)
		elif len(nod.loc)==3:
			xyz.InsertPoint(nod.localID, nod.loc[0]+nodDof[0],  nod.loc[1]+nodDof[1], nod.loc[2]+nodDof[2])
			uvw.SetTuple3(nod.localID, nodDof[0], nodDof[1], nodDof[2])
		else:
			raise "node.loc has wrong length for plotting"

		if len(nodDof)>len(nod.loc):
			rho.SetValue(nod.localID, nodDof[-1]);
		else:
			rho.SetValue(nod.localID,-1);

	# Initialize grid and element types
	ug = vtk.vtkUnstructuredGrid()
	ug.SetPoints(xyz)
	
	line_type = vtk.vtkLine().GetCellType()
	quad_type = vtk.vtkQuad().GetCellType()
	hexa_type = vtk.vtkHexahedron().GetCellType()
	
	# Extract element connectivities and element output variables
	count = 0;
	idList = vtk.vtkIdList(); idList.Allocate(8)
	for el in mc.elements:
		if(el.type=='CBAR' or el.type=='CBAR1' or el.type=='CBARX' or el.type=='CBEAM' or el.type=='CBEAMX'):
			el_type = line_type
			Jg.SetValue(count, 1.)
			Stretch.SetValue(count,eh.getStretch(el,0,dc))
			
		elif(el.type=='CQUAD' or el.type=='CQUAD7'):
			el_type = quad_type
			
			Jg.SetValue(count, np.mean([eh.getJg(el,ip,dc,2) for ip in range(4)]))
			Stretch.SetValue(count,1.)
		elif(el.type=='CHEXA'or el.type=='CHEXA7'):
			el_type = hexa_type

			Jg.SetValue(count,np.mean([eh.getJg(el,ip,dc,3) for ip in range(8)]))
			Stretch.SetValue(count,1.)
		else:
			raise "Element type "+el.type+" is not supported for output plotting"


		idList.Reset()
		[idList.InsertId(i,el.nodes[i].localID) for i in range(len(el.nodes))]
		ug.InsertNextCell(el_type, idList)

		PropID.SetValue(count, el.property.localID)

		try:
			State.SetValue(count, el.state)
		except :
			State.SetValue(count, -1)
		count = count+1

	# Add data
	ug.GetPointData().AddArray(uvw)
	ug.GetPointData().AddArray(rho)
	ug.GetCellData().AddArray(Jg)
	ug.GetCellData().AddArray(Stretch)
	ug.GetCellData().AddArray(PropID)
	ug.GetCellData().AddArray(State)
	
	# Write to file
	#filename = 'paraview_'+'{0:04d}'.format(dc.step)+'{0:04d}'.format(it)+'.vtu'
	filename = 'paraview_'+'{0:04d}'.format(dc.step)+'.vtu'
	w = vtk.vtkXMLUnstructuredGridWriter()
	#w.SetInput(ug)
	w.SetInputData(ug)
	w.SetFileName(mc.folderPostProcess+'/'+filename)
	w.Write()

	# Create paraview.pvd file if not existent
	cdef object pvd_file
	pvdFilename = mc.folderPostProcess+'/paraview.pvd'
	if not os.path.exists(pvdFilename):
		pvd_file=open(pvdFilename, 'w')
		pvd_file.write("<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">\n\t<Collection>\n\t\t#REPLACE_HERE#\n\t</Collection>\n</VTKFile>")
		pvd_file.close()
	
	# Update paraview.pvd file
	pvd_file = open(pvdFilename,'r')
	fileData = pvd_file.read()
	pvd_file.close()

	pvd_file = open(pvdFilename,'w')
	pvd_file.write(fileData.replace("#REPLACE_HERE#","<DataSet part=\"0\" timestep=\""+str(dc.time)+"\" file=\""+filename+"\"/>\n\t\t#REPLACE_HERE#"))
	pvd_file.close()



cpdef PlotOutput(ModelContainer mc, DataContainer dc):
	cdef Element el
	cdef list img
	cdef ndim = mc.numDofPerNode
	img = []
	
	fig, ax = plt.subplots()
	for el in mc.elements:
		if(len(el.nodes)==4):
			PlotFace(dc,el,el.nodes,ndim)
		elif(len(el.nodes)==8):
			PlotFace(dc,el,[el.nodes[0],el.nodes[1],el.nodes[2],el.nodes[3]],ndim) # bottom face
			PlotFace(dc,el,[el.nodes[4],el.nodes[5],el.nodes[6],el.nodes[7]],ndim) # top face
			PlotFace(dc,el,[el.nodes[0],el.nodes[1],el.nodes[5],el.nodes[4]],ndim) # front face
			PlotFace(dc,el,[el.nodes[2],el.nodes[3],el.nodes[7],el.nodes[6]],ndim) # rear face
			PlotFace(dc,el,[el.nodes[0],el.nodes[3],el.nodes[7],el.nodes[4]],ndim) # left face
			PlotFace(dc,el,[el.nodes[1],el.nodes[2],el.nodes[6],el.nodes[5]],ndim) # right face

	ax.text(0.0, 0.1, 'Time = %.2f [days]'%(dc.time), transform=ax.transAxes, fontsize=14,
        verticalalignment='bottom')

	plt.axis('equal')
	plt.axis('off')

	plt.savefig(mc.folderPostProcess+'/eps_'+'{0:04d}'.format(dc.step)+'.eps',format='eps')
	plt.savefig(mc.folderPostProcess+'/png_'+'{0:04d}'.format(dc.step)+'.png',dpi=100,format='png')


	plt.close()

cpdef PlotFace(DataContainer dc, Element el, list n, int ndim):
	cdef np.ndarray[double, ndim=1] faceColor
	cdef np.ndarray[double, ndim=2] ndof

	# Get dof at nodes
	ndof = np.array([n[0].Dof(dc), n[1].Dof(dc), n[2].Dof(dc), n[3].Dof(dc)])

	# Compute faceColor based on property
	if(el.property.localID<200):						# Gray matter
		faceColor = np.array([255., 103., 103.])/255.
	elif(el.property.localID<300):						# White matter
		faceColor = np.array([255., 204., 153.])/255.
	else:												# Intermediate layer
		faceColor = np.array([93.,  74.,  56.])/255. 

	if(ndim==2):
		plt.fill([n[0].x+ndof[0,0], n[1].x+ndof[1,0], n[2].x+ndof[2,0], n[3].x+ndof[3,0]],\
				 [n[0].y+ndof[0,1], n[1].y+ndof[1,1], n[2].y+ndof[2,1], n[3].y+ndof[3,1]],\
				 fc=faceColor,ec=[0.5,0.5,0.5],lw=0.2)
	elif(ndim==3):
		plt.fill([n[0].x+ndof[0,0], n[1].x+ndof[1,0], n[2].x+ndof[2,0], n[3].x+ndof[3,0]],\
				 [n[0].y+ndof[0,1], n[1].y+ndof[1,1], n[2].y+ndof[2,1], n[3].y+ndof[3,1]],\
				 fc=faceColor,ec=[0.5,0.5,0.5],lw=0.2)

cpdef WriteToLog(ModelContainer mc, str text):
	cdef object log_file
	print text
	log_file=open(mc.folderOutput + '/log.txt', "a")
	log_file.write(text+"\n")
	log_file.close()

cpdef WriteToOutput(ModelContainer mc, str filename, str text):
	cdef object out_file
	out_file=open(mc.folderOutput + '/'+filename, "a")
	out_file.write(text+"\n")
	out_file.close()
