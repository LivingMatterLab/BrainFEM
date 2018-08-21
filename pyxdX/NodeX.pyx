# -*- coding: utf-8 -*-
from Node cimport *
cimport numpy as np 
import numpy as np

cdef class NodeX(Node):
    def __init__(self,loc):
        super(NodeX,self).__init__(loc);

        # Initialize element towards plus and minus side
        self.elMinus = None
        self.elPlus = None

        # Initialize mpc that governs this node
        self.mpc = None

        # Initialize curvature at this node
        self.curvature0 = 0.0;

    def __str__(self):
        cdef str eidM, eidP, mid

        eidM = "None" if self.elMinus==None else str(self.elMinus.localID)
        eidP = "None" if self.elPlus==None else str(self.elPlus.localID)
        mid = "None" if self.mpc==None else str(self.mpc.slaveDofID)

        return super(NodeX,self).__str__() + "\t [el-,el+] = [" + eidM  \
                + "," + eidP+"] \t mpcSlaveDofID = " + mid

    cpdef InitializeCurvature(self,DataContainer dc):
        cdef np.ndarray[double,ndim=1] locM, loc, locP,kappa
        cdef double lM, lP, gamma

        # Location of node on minus side (locM), this node (loc), and node on plus side (locP)
        # 'if statements' account for far left and far right node
        if(self.elMinus==None):
            locM = np.array(self.loc)                  + np.array(self.Dof(dc))
            loc  = np.array(self.elPlus.nodes[1].loc)  + np.array(self.elPlus.nodes[1].Dof(dc))
            #locP = np.array(self.elPlus.nodes[1].elPlus.nodes[1].loc) + np.array(self.elPlus.nodes[1].elPlus.nodes[1].Dof(dc))    
            locP = 2.*loc-locM  # Ensure curvature=0
        elif(self.elPlus==None):
            locM = np.array(self.elMinus.nodes[0].elMinus.nodes[0].loc) + np.array(self.elMinus.nodes[0].elMinus.nodes[0].Dof(dc))
            loc = np.array(self.elMinus.nodes[0].loc) + np.array(self.elMinus.nodes[0].Dof(dc))
            #locP  = np.array(self.loc)                  + np.array(self.Dof(dc))
            locP = 2.*loc-locM  # Ensure curvature=0
        else:
            locM = np.array(self.elMinus.nodes[0].loc) + np.array(self.elMinus.nodes[0].Dof(dc))
            loc  = np.array(self.loc)                  + np.array(self.Dof(dc))
            locP = np.array(self.elPlus.nodes[1].loc)  + np.array(self.elPlus.nodes[1].Dof(dc))

        lM = np.linalg.norm(loc-locM)           # Current length of element on the minus side
        lP = np.linalg.norm(locP-loc)           # Current length of element on the plus side
        gamma = 1./(0.5*lM*lP*(lM+lP))          # Parameter to compute curvature

        kappa = gamma*(lP*locM-(lP+lM)*loc+lM*locP) # Curvature vector at this node


        self.curvature0 = np.linalg.norm(kappa)      # Norm of kappa, i.e. curvature scalar

        """
        # Direction of curvature
        if(self.curvature<1e-6):
            self.curv_direction = np.zeros(len(self.loc))
        else:
            self.curv_direction = kappa/self.curvature0
        """

        """
        print "Computing curvature of node ", self.localID
        
        print "locM  = " , locM
        print "loc   = " , loc
        print "locP  = " , locP
        print "lM    = " , lM
        print "lP    = " , lP
        print "gamma = " , gamma
        print "kappa = " , kappa
        
        print "curvature0      = ", self.curvature0
        """
        