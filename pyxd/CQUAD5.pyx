# -*- coding: utf-8 -*-
# CQUAD5 - axisymmetric around z=0 axis (coordinates are in R,Z plane)
from CQUAD cimport *
from math import *
cimport ElementHelper as eh

cimport numpy as np
import numpy as np
import scipy.linalg as spl

cdef class CQUAD5(CQUAD):
    def __init__(self,nodes,elementDirection,elementProperty):
        super(CQUAD5,self).__init__(nodes,elementDirection,elementProperty);

        # Compute R0 at the integration points
        self.ComputeR0()


    def __str__(self):
        return super(CQUAD5,self).__str__()
    
    cpdef InitializeData(self,DataContainer dc):
        # Initialize Fg and Fg0
        fgLoc = [];
        for i in range(0,4):
            fgLoc.append(np.eye(3))

        self.Fg0ID =  dc.AddToData(dc.Fg0,dc.countFg0,fgLoc)
        self.FgID  =  dc.AddToData(dc.Fg,dc.countFg,fgLoc)
    
    cpdef ComputeR0(self):
        cdef int ip, ndim, i
        cdef double r_ip, s_ip
        cdef np.ndarray[double, ndim=1] r0,N
        cdef np.ndarray[double, ndim=2] gp
        
        r0  = self.Loc();
        gp = 0.577350269189626* \
             np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])
        ndim = 2;

        # Loop through integration points
        self.R0 = np.zeros(4)
        for ip in range(4):
            # r,s coordinates at this integration point
            r_ip = gp[ip,0];
            s_ip = gp[ip,1];

            N = eh.Jacobian_Q2(r_ip,s_ip,r0)['N']

            for i in range(4):
                self.R0[ip]+=N[i]*r0[ndim*i]

    cpdef BuildElementMatrix(self,DataContainer dc):
        cdef int ndim,ndimF, ip, i,j,p,q,r,s
        cdef double rc,r_ip,s_ip, detJ
        cdef np.ndarray[double,ndim=1] ve, veE, r0, r0E, wp,N
        cdef np.ndarray[double,ndim=2] gp, mat, dNx, J,F
        cdef np.ndarray[double,ndim=3] B
        cdef np.ndarray[double,ndim=4] A
        cdef dict deform, M

        ve  = self.Dof(dc);
        veE = self.T.dot(ve);        # ve in element coordinates
        r0  = self.Loc();
        r0E = self.T.dot(r0);        # r0 in element coordinates

        gp = 0.577350269189626* \
             np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])

        wp = np.array([1.,1.,1.,1.]);
        mat = np.zeros((len(ve),len(ve)));
        ndim = 2;
        ndimF = 3;

        # Loop through integration points
        for ip in range(4):
            # r,s coordinates at this integration point
            r_ip = gp[ip,0];
            s_ip = gp[ip,1];

             # Compute deformation gradient
            deform = eh.LocalDeformation_Q2(r_ip,s_ip,r0E,veE)
            dNx = deform['dNx'].T
            N   = deform['N']
            J = deform['J'];
            detJ = np.linalg.det(J)

            # Compute current radius
            rc = 0.
            for i in range(4):
                rc+=N[i]*(r0[ndim*i]+ve[ndim*i])

            # Compute tangent matrix and updated deformation gradient
            F = spl.block_diag(deform['F'],np.array([rc/self.R0[ip]]))
            M = self.property.Piola1Stiffness(F,eh.getFg0(self,ip,dc,ndimF),dc.dt,ndimF);
            A = M['A'];
            eh.setFg(self,ip,M['Fg'],dc,ndimF)
            
            for i in range(4):
                for j in range(i+1):
                    eni = i*ndim;
                    enj = j*ndim;

                    for p in range(ndim):
                        for q in range(ndim):
                            for r in range(ndim):
                                for s in range(ndim):
                                    mat[eni+p,enj+q]+= \
                                                    2.*pi*dNx[r,i]*A[p,r,q,s]*dNx[s,j]*detJ*wp[ip]*self.R0[ip];

                                mat[eni+p,enj+q]+= 2.*pi*N[i]*(p==0)*A[ndim,ndim,q,r]*dNx[r,j]/self.R0[ip]*detJ*wp[ip]*self.R0[ip];
                                mat[eni+p,enj+q]+= 2.*pi*N[j]*(q==0)*A[p,r,ndim,ndim]*dNx[r,i]/self.R0[ip]*detJ*wp[ip]*self.R0[ip];
                            
                            mat[eni+p,enj+q]+= 2.*pi*N[i]*(p==0)*A[ndim,ndim,ndim,ndim]*N[j]*(q==0)/self.R0[ip]**2*detJ*wp[ip]*self.R0[ip];
                      
        
        # Make symmetric
        mat = np.tril(mat)+np.tril(mat,-1).T;

        # Transform back into global x,y coordinates
        mat = self.T.T.dot(mat).dot(self.T);

        return mat

    cpdef BuildInternalForceVector(self,DataContainer dc):
        cdef int ndim,ndimF, ip, i,p,q
        cdef double rc,r_ip,s_ip, detJ
        cdef np.ndarray[double,ndim=1] ve, veE, r0, r0E,N, wp, vec
        cdef np.ndarray[double,ndim=2] gp, dNx, J, P, F
        cdef dict deform, M

        ve  = self.Dof(dc);
        veE = self.T.dot(ve);        # ve in element coordinates
        r0  = self.Loc();
        r0E = self.T.dot(r0);        # r0 in element coordinates

        gp = 0.577350269189626* \
             np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])

        wp = np.array([1.,1.,1.,1.]);
        vec = np.zeros(len(ve));
        ndim = 2;
        ndimF = 3;

        # Loop through integration points
        for ip in range(4):
            # r,s coordinates at this integration point
            r_ip = gp[ip,0];
            s_ip = gp[ip,1];

            # Compute deformation gradient
            deform = eh.LocalDeformation_Q2(r_ip,s_ip,r0E,veE)
            dNx = deform['dNx'].T
            N   = deform['N']
            J = deform['J'];
            detJ = np.linalg.det(J)

            # Compute current radius
            rc = 0.
            for i in range(4):
                rc+=N[i]*(r0[ndim*i]+ve[ndim*i])

            # Compute tangent matrix and updated deformation gradient
            F = spl.block_diag(deform['F'],np.array([rc/self.R0[ip]]))
            M = self.property.Piola1Stiffness(F,eh.getFg(self,ip,dc,ndimF),0,ndimF);
            P = M['P'];

            # Check F
            Fcheck = np.zeros((ndimF,ndimF))
            
            for i in range(4):
                eni = i*ndim;

                for p in range(ndim):
                    for q in range(ndim):
                        vec[eni+p] += 2*pi*P[p,q]*dNx[q,i]*detJ*wp[ip]*self.R0[ip];
                    vec[eni+p]+= 2*pi*N[i]*P[ndim,ndim]*(p==0)/self.R0[ip]*detJ*wp[ip]*self.R0[ip];
                
    
        # Transform back into global x,y,z coordinates
        vec = self.T.T.dot(vec)

        return vec

