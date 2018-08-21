from Solver cimport *
cdef class ModelContainer:
	cdef public list materials, properties, amplitudes, nodes, elements, spc, mpc,loads,contacts,pressures,mechanisms
	cdef public str folderMain, folderSimulation, folderInput, folderOutput, folderPostProcess
	cdef public int nmat, nprop, nampl, nnode, nel, nspc, nmpc, nloads, ncontacts,npressures, nmech
	cdef public int numDispDofPerNode, numDofPerNode, ndof
	cpdef SetListLengths(self)
	cpdef BuildModel(self,object p)
	cpdef CreateOutputFolders(self)
	cpdef WriteStepOutput(self, Solver s)



