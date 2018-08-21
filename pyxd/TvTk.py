from tvtk.api import tvtk


#Generate Paraview file
cpdef ParaviewOutput(ModelContainer mc, DataContainer dc):
	cdef np.ndarray[double,ndim=1] Jg, PropID
	cdef np.ndarray[double,ndim=2] xyz, uvw
	cdef np.ndarray[np.int_t,ndim=1] State

	xyz=np.zeros((mc.nnode,3))
	uvw=np.zeros((mc.nnode,3))
	Jg     = np.zeros(mc.nel) 						# Growth factor
	State  = np.zeros(mc.nel).astype(np.int)		# Element state
	PropID = np.zeros(mc.nel) 						# Property ID

	for nod in mc.nodes:
		if len(nod.loc)==2:	
			xyz[nod.localID]= np.array([nod.loc[0]+nod.Dof(dc)[0],nod.loc[1]+nod.Dof(dc)[1],0.])
			uvw[nod.localID]= np.array([nod.Dof(dc)[0],nod.Dof(dc)[1],0.])
		elif len(nod.loc)==3:
			xyz[nod.localID]= np.array([nod.loc[0]+nod.Dof(dc)[0],nod.loc[1]+nod.Dof(dc)[1],nod.loc[2]+nod.Dof(dc)[2]])
			uvw[nod.localID]= np.array([nod.Dof(dc)[0],nod.Dof(dc)[1],nod.Dof(dc)[2]])
		else:
			raise "node.loc has wrong length for plotting"

	ug = tvtk.UnstructuredGrid(points=xyz)
	line_type = tvtk.Line().cell_type
	quad_type = tvtk.Quad().cell_type
	hexa_type = tvtk.Hexahedron().cell_type
	
	count = 0;
	for el in mc.elements:
		if(el.type=='CBAR' or el.type=='CBAR1' or el.type=='CBARX' or el.type=='CBEAMX'):
			el_type = line_type
			Jg[count]= 1.
			
		elif(el.type=='CQUAD'):
			el_type = quad_type
			
			Jg[count]=np.mean([eh.getJg(el,ip,dc,2) for ip in range(4)])
		elif(el.type=='CHEXA'):
			el_type = hexa_type

			Jg[count]=np.mean([eh.getJg(el,ip,dc,3) for ip in range(8)])
		else:
			raise "Element type "+el.type+" is not supported for output plotting"

		ug.insert_next_cell(el_type,np.array([nod.localID for nod in el.nodes]).flatten())
		PropID[count]= el.property.localID

		try:
			State[count] = el.state;
		except :
			State[count] = -1;
		count = count+1

	
	pd        = tvtk.PointData()
	"""
	pd.scalars = v
	pd.scalars.name = 'v'
	ug.point_data.add_array(pd.scalars)
	"""

	pd.vectors = uvw.tolist()
	pd.vectors.name = 'Disp'
	ug.point_data.add_array(pd.vectors)

	cdJg        = tvtk.CellData()
	cdJg.scalars = Jg
	cdJg.scalars.name='J_g'
	ug.cell_data.add_array(cdJg.scalars)
	

	cdPID        = tvtk.CellData()
	cdPID.scalars = PropID
	cdPID.scalars.name='PropID'
	ug.cell_data.add_array(cdPID.scalars)

	cdState        = tvtk.CellData()
	cdState.scalars = State
	cdState.scalars.name='State'
	ug.cell_data.add_array(cdState.scalars)

	filename = mc.folderPostProcess+'/paraview_'+'{0:04d}'.format(dc.step)+'.vtu'
	w = tvtk.XMLUnstructuredGridWriter(input = ug, file_name=(filename))
	w.write()