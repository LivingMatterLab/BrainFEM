# -*- coding: utf-8 -*-
cimport numpy as np
import numpy as np
cdef class Element(object):
    def __init__(self,elementType,nodes,elementDirection,elementProperty):
        self.localID = -1;
        self.type = elementType;
        self.nodes = nodes;
        self.direction = elementDirection;
        self.property = elementProperty;

    cpdef Loc(self):
        cdef list loc
        loc = [];
        for nod in self.nodes:
            loc.extend(nod.loc)
        return np.array(loc)

    cpdef DofID(self):
        cdef list dofID
        dofID = [];
        for nod in self.nodes:
            dofID.extend(nod.dofID)
        return np.array(dofID)

    cpdef Dof0(self,DataContainer dc):
        return np.array([dc.dof0[i] for i in self.DofID()])
        
    cpdef Dof(self,DataContainer dc):
        return np.array([dc.dof[i] for i in self.DofID()])


    def __str__(self):
        return "Element " + str(self.localID)+":\t type = " + self.type +"\t prop = " + str(self.property.localID) +",\t nodes = [" + ', '.join(str(nod.localID) for nod in self.nodes) + "]"
