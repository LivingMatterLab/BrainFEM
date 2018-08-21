# -*- coding: utf-8 -*-
cimport numpy as np
import numpy as np

cdef class Node:
    def __init__(self, loc):
        self.loc = loc
        self.x = loc[0]
        self.y = loc[1]
        try:
            self.z = loc[2]
        except:
            self.z = 0.
        self.dofID = [];        # ID of nodal dof in global matrices/vectors      
        self.localID = -1;
        
    def __str__(self):
        return "Node " + str(self.localID)+":  \t loc = ["+ ', '.join(str(x) for x in self.loc)+"]" \
                       +"  \t dofID = ["+ ', '.join(str(x) for x in self.dofID)+"]"

    cpdef Dof(self,DataContainer dc):
        return [dc.dof[i] for i in self.dofID]
    
    cpdef Dof0(self,DataContainer dc):
        return [dc.dof0[i] for i in self.dofID]

    cpdef DispDof(self,DataContainer dc):
        return [dc.dof[self.dofID[i]] for i in range(min(len(self.loc),len(self.dofID)))]

    cpdef DispDof0(self,DataContainer dc):
        return [dc.dof0[self.dofID[i]] for i in range(min(len(self.loc),len(self.dofID)))]