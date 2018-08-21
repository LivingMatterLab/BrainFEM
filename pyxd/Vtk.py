
import vtk
cpdef ParaviewOutput(ModelContainer mc, DataContainer dc):
	cdef list nodDof

	# Initialize output variables
	Jg = vtk.vtkDoubleArray()
	Jg.SetName("Jg")
	Jg.SetNumberOfTuples(mc.nel)

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
		if(el.type=='CBAR' or el.type=='CBAR1' or el.type=='CBARX' or el.type=='CBEAMX'):
			el_type = line_type
			Jg.SetValue(count, 1.)
			
		elif(el.type=='CQUAD'):
			el_type = quad_type
			
			Jg.SetValue(count, np.mean([eh.getJg(el,ip,dc,2) for ip in range(4)]))
		elif(el.type=='CHEXA'):
			el_type = hexa_type

			Jg.SetValue(count,np.mean([eh.getJg(el,ip,dc,3) for ip in range(8)]))
		else:
			raise "Element type "+el.type+" is not supported for output plotting"

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
	ug.GetCellData().AddArray(Jg)
	ug.GetCellData().AddArray(PropID)
	ug.GetCellData().AddArray(State)
	
	# Write to file
	filename = mc.folderPostProcess+'/paraview_'+'{0:04d}'.format(dc.step)+'.vtu'
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInput(ug)
	#w.SetInputData(ug)
	w.SetFileName(filename)
	w.Write()