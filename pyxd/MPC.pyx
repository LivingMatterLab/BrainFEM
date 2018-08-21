# -*- coding: utf-8 -*-
cdef class MPC:
    def __init__(self, nodes, dofIDs,weights):
        self.slaveDofID     = nodes[0].dofID[dofIDs[0]]
        self.masterDofIDs   = [nodes[i].dofID[dofIDs[i]] for i in range(1,len(nodes))]
        self.masterWeights  = [weights[i] for i in range(1,len(nodes))]
        
        
    def __str__(self):
        return "MPC: \t SlaveDofID " + str(self.slaveDofID)+"  \t MasterDofIDs = ["+ ', '.join(str(x) for x in self.masterDofIDs)+"]" \
                       +"\t MasterWeights = ["+ ', '.join(str(x) for x in self.masterWeights)+"]"
