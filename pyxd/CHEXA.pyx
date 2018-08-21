# -*- coding: utf-8 -*-
from Element cimport *
from math import *
cimport ElementHelper as eh

cimport numpy as np
import numpy as np
import scipy.linalg as spl

cdef class CHEXA(Element):
    def __init__(self,nodes,elementDirection,elementProperty):
        super(CHEXA,self).__init__('CHEXA',nodes,elementDirection,elementProperty);

        # Compute element rotation matrix T
        self.T = self.RotationMatrix()

    def __str__(self):
        return super(CHEXA,self).__str__()

    cpdef InitializeData(self,DataContainer dc):
        # Initialize Fg and Fg0
        fgLoc = [];
        for i in range(0,8):
            fgLoc.append(np.eye(3))

        self.Fg0ID =  dc.AddToData(dc.Fg0,dc.countFg0,fgLoc)
        self.FgID  =  dc.AddToData(dc.Fg,dc.countFg,fgLoc)     

    cpdef RotationMatrix(self):
        cdef np.ndarray[double, ndim=1] ex,ez, ey
        cdef np.ndarray[double, ndim=2] nodeT

        ex = np.array(self.direction);
        ez = np.array([0.,0.,1.]);
        ez = ez - np.dot(ez,ex)*ex;
        if(np.linalg.norm(ez)<0.01):
            ez = np.array([1.,0.,0.]);
            ez = ez - np.dot(ez,ex)*ex;
        ez = ez/np.linalg.norm(ez);
        ey = np.cross(ez,ex);

        nodeT = np.squeeze([ex,ey,ez])
        return spl.block_diag(nodeT,nodeT,nodeT,nodeT,nodeT,nodeT,nodeT,nodeT)

    def BuildElementMatrixPublic(self,DataContainer dc):
        return self.BuildElementMatrix(dc)

    cpdef BuildElementMatrix(self,DataContainer dc):
        cdef int ndim, ip, i,j,p,q,r,s
        cdef double r_ip,s_ip,t_ip, detJ
        cdef np.ndarray[double,ndim=1] ve, veE, r0, r0E, wp
        cdef np.ndarray[double,ndim=2] gp, mat, dNx, J
        cdef np.ndarray[double,ndim=4] A
        cdef dict deform, M

        ve  = self.Dof(dc);
        veE = self.T.dot(ve);        # ve in element coordinates
        r0  = self.Loc();
        r0E = self.T.dot(r0);        # r0 in element coordinates

        gp = 0.577350269189626* \
             np.array([[-1., -1., -1.],[1., -1., -1.],[1., 1., -1.], [-1., 1., -1.],\
                       [-1., -1., 1.], [1., -1., 1.],[1., 1., 1.],[-1., 1., 1.]])

        wp = np.array([1.,1.,1.,1.,1.,1.,1.,1.]);
        mat = np.zeros((len(ve),len(ve)));
        ndim = 3;

        # Loop through integration points
        for ip in range(8):
            # r,s,t coordinates at this integration point
            r_ip = gp[ip,0];
            s_ip = gp[ip,1];
            t_ip = gp[ip,2];

            # Compute deformation gradient
            deform = eh.LocalDeformation_H3(r_ip,s_ip,t_ip,r0E,veE)

            # Compute tangent matrix and updated deformation gradient
            M = self.property.Piola1Stiffness(deform['F'],eh.getFg0(self,ip,dc,ndim),dc.dt,ndim);
            A = M['A'];
            eh.setFg(self,ip,M['Fg'],dc,ndim)

            dNx = deform['dNx'].T
            J = deform['J'];
            detJ = np.linalg.det(J)
            
            for i in range(8):
                for j in range(i+1):
                    eni = i*ndim;
                    enj = j*ndim;

                    for p in range(ndim):
                        for q in range(ndim):
                            for r in range(ndim):
                                for s in range(ndim):
                                    mat[eni+p,enj+q] = mat[eni+p,enj+q] \
                                                       + dNx[r,i]*A[p,r,q,s]*dNx[s,j]*detJ*wp[ip];

        
        # Make symmetric
        mat = np.tril(mat)+np.tril(mat,-1).T;

        # Transform back into global x,y,z coordinates
        mat = self.T.T.dot(mat).dot(self.T);

        return mat

    def BuildInternalForceVectorPublic(self,DataContainer dc):
        return self.BuildInternalForceVector(dc)

    cpdef BuildInternalForceVector(self,DataContainer dc):
        cdef int ndim, ip, i,p,q
        cdef double r_ip,s_ip,t_ip, detJ
        cdef np.ndarray[double,ndim=1] ve, veE, r0, r0E, wp, vec
        cdef np.ndarray[double,ndim=2] gp, dNx, J, P
        cdef dict deform, M

        ve  = self.Dof(dc);
        veE = self.T.dot(ve);        # ve in element coordinates
        r0  = self.Loc();
        r0E = self.T.dot(r0);        # r0 in element coordinates

        gp = 0.577350269189626* \
             np.array([[-1., -1., -1.],[1., -1., -1.],[1., 1., -1.], [-1., 1., -1.],\
                       [-1., -1., 1.], [1., -1., 1.],[1., 1., 1.],[-1., 1., 1.]])

        wp = np.array([1.,1.,1.,1.,1.,1.,1.,1.]);
        vec = np.zeros(len(ve));
        ndim = 3;

        # Loop through integration points
        for ip in range(8):
            # r,s,t coordinates at this integration point
            r_ip = gp[ip,0];
            s_ip = gp[ip,1];
            t_ip = gp[ip,2];

            # Compute deformation gradient
            deform = eh.LocalDeformation_H3(r_ip,s_ip,t_ip,r0E,veE)

            # Compute tangent matrix and updated deformation gradient
            M = self.property.Piola1Stiffness(deform['F'],eh.getFg(self,ip,dc,ndim),0,ndim);
            P = M['P'];

            dNx = deform['dNx'].T
            J = deform['J'];
            detJ = np.linalg.det(J)
            
            for i in range(8):
                eni = i*ndim;
                for p in range(ndim):
                    for q in range(ndim):
                        vec[eni+p] = vec[eni+p] \
                                           + P[p,q]*dNx[q,i]*detJ*wp[ip];

        # Transform back into global x,y,z coordinates
        vec = self.T.T.dot(vec)

        return vec