# -*- coding: utf-8 -*-
# CQUAD requires that the problem is 2D, i.e. only 2 dof per node
from Element cimport *
from math import *
cimport ElementHelper as eh

cimport numpy as np
import numpy as np
import scipy.linalg as spl

cdef class CQUAD7(Element):
    def __init__(self,nodes,elementDirection,elementProperty):
        super(CQUAD7,self).__init__('CQUAD7',nodes,elementDirection,elementProperty);

        # Compute element rotation matrix T
        self.T = self.RotationMatrix()

    def __str__(self):
        return super(CQUAD7,self).__str__()

    cpdef InitializeData(self,DataContainer dc):
        # Initialize Fg, Fg0, nc, nc0
        fgLoc = [];
        ncLoc = [];
        for i in range(0,4):
            fgLoc.append(np.eye(2))
            ncLoc.append(0.);

        self.Fg0ID =  dc.AddToData(dc.Fg0,dc.countFg0,fgLoc)
        self.FgID  =  dc.AddToData(dc.Fg,dc.countFg,fgLoc)

        self.nc0ID =  dc.AddToData(dc.nc0,dc.countNc0,ncLoc)
        self.ncID  =  dc.AddToData(dc.nc,dc.countNc,ncLoc)

    cpdef RotationMatrix(self):
        cdef double th
        cdef np.ndarray[double, ndim=1] ex, ey
        cdef np.ndarray[double, ndim=2] nodeT

        ex = self.direction;
        
        ex = np.array(self.direction)
        th = np.arctan2(ex[1],ex[0]);
        ey = np.array([-sin(th), cos(th)]);

        nodeT = np.squeeze([ex,ey])
        return spl.block_diag(nodeT,np.diag([1]),nodeT,np.diag([1]),\
                              nodeT,np.diag([1]),nodeT,np.diag([1]))

    cpdef BuildElementMatrix(self,DataContainer dc):
        cdef int ndim, ip, i,j,p,q,r,s
        cdef double r_ip,s_ip, detJ, Rho
        cdef double detF, detF0, detF_dt
        cdef np.ndarray[np.int64_t,ndim=1] udof, rdof
        cdef np.ndarray[double,ndim=1] ve, veE, ve0,ve0E, r0, r0E,N, wp,gradRho, Q
        cdef np.ndarray[double,ndim=2] gp, mat, dNx, J, Finv, d, dfdF, dPdRho
        cdef np.ndarray[double,ndim=3] D3
        cdef np.ndarray[double,ndim=4] A
        cdef dict deform, M, Diff, RhoS

        udof = np.array([0, 1, 3, 4, 6, 7, 9, 10]);
        rdof = np.array([2, 5, 8, 11]);

        ve  = self.Dof(dc);
        veE = self.T.dot(ve);        # ve in element coordinates
        ve0  = self.Dof0(dc);
        ve0E = self.T.dot(ve0);        # ve0 in element coordinates
        r0  = self.Loc();
        r0E = self.T[np.ix_(udof,udof)].dot(r0);        # r0 in element coordinates

        gp = 0.577350269189626* \
             np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])

        wp = np.array([1.,1.,1.,1.]);
        mat = np.zeros((len(ve),len(ve)));
        ndim = 2;

        # Loop through integration points
        for ip in range(4):
            # r,s coordinates at this integration point
            r_ip = gp[ip,0];
            s_ip = gp[ip,1];

            # Compute deformation gradient
            deform = eh.LocalDeformation_Q2(r_ip,s_ip,r0E,veE[udof])
            N = deform['N']
            dNx = deform['dNx'].T
            J = deform['J'];
            detJ = np.linalg.det(J)
            Finv = np.linalg.inv(deform['F']);

            # Compute deformation gradient
            deform0 = eh.LocalDeformation_Q2(r_ip,s_ip,r0E,ve0E[udof])
            detF0 = np.linalg.det(deform0['F']);
            detF  = np.linalg.det(deform['F']);
            detF_dt = (detF-detF0)/dc.dt

            # Compute gradient of rho wrt reference X
            Rho = 0;
            gradRho = np.zeros(ndim);
            for j in range(4):
                Rho = Rho + N[j]*veE[rdof[j]];
                for s in range(ndim):
                     gradRho[s] += dNx[s,j]*veE[rdof[j]];
            # Compute gradient of rho wrt current x
            gradRho = Finv.T.dot(gradRho);

            # Compute tangent matrix and updated deformation gradient
            M = self.property.Piola1Stiffness(deform['F'],eh.getFg0(self,ip,dc,ndim), Rho, eh.getNc0(self,ip,dc),dc.dt,ndim);
            Diff = self.property.Diffusivity(deform['F'],Rho,gradRho,ndim);
            RhoS = self.property.SourceRho(deform['F'],Rho,gradRho,dc,ndim);
            A = M['A'];
            dPdRho = M['dPdRho'];
            Q = Diff['q']
            d = Finv.dot(Diff['D2'].dot(Finv.T));
            d1 = Finv.dot(Diff['D1'])
            D3 = Diff['D3']
            dfdF = RhoS['dfdF'];
            eh.setFg(self,ip,M['Fg'],dc,ndim)
            eh.setNc(self,ip,M['nc'],dc)
            
            for i in range(4):
                for j in range(4):
                    eni = i*(ndim+1);        # +1 due to additional field
                    enj = j*(ndim+1);        # +1 due to additional field

                    # Displacement contribution
                    for p in range(ndim):
                        for q in range(ndim):
                            for r in range(ndim):
                                for s in range(ndim):
                                    mat[eni+p,enj+q] += dNx[r,i]*A[p,r,q,s]*dNx[s,j]*detJ*wp[ip];

                    # Density contribution
                    mat[eni+ndim,enj+ndim] += N[i]*N[j]*(1./dc.dt+detF_dt/detF)*detJ*wp[ip]
                    mat[eni+ndim,enj+ndim] -= N[i]*N[j]*RhoS['dfdRho']*detJ*wp[ip]
                    for r in range(ndim):
                        mat[eni+ndim,enj+ndim]+= dNx[r,i]*d1[r]*N[j]*detJ*wp[ip];
                        for s in range(ndim):
                            mat[eni+ndim,enj+ndim]+= dNx[r,i]*d[r,s]*dNx[s,j]*detJ*wp[ip];
                    
                    # coupling contribution (dRrho/dPhi)
                    for p in range(ndim):
                        for q in range(ndim):
                            mat[eni+ndim,enj+p] += Rho*detF0/detF/dc.dt*N[i]*Finv[q,p]*dNx[q,j]*detJ*wp[ip]
                            mat[eni+ndim,enj+p] -= N[i]*dfdF[p,q]*dNx[q,j]*detJ*wp[ip]
                            for r in range(ndim):
                                for s in range(ndim):
                                    mat[eni+ndim,enj+p] += dNx[r,i]*(-Finv[r,p]*Finv[s,q]*Q[q]+Finv[r,q]*D3[q,p,s])*dNx[s,j]*detJ*wp[ip]

                    # coupling contribution (dRphi/dRho)
                    for p in range(ndim):
                        for r in range(ndim):
                            mat[eni+p,enj+ndim] += dNx[r,i]*dPdRho[p,r]*N[j]*detJ*wp[ip]
                            
        # Make symmetric
        #mat = np.tril(mat)+np.tril(mat,-1).T;

        # Transform back into global x,y coordinates
        mat = self.T.T.dot(mat).dot(self.T);
        return mat

    cpdef BuildInternalForceVector(self,DataContainer dc):
        cdef int ndim, ip, i,p,q,s
        cdef double r_ip,s_ip, detJ, Rho,dRhodt
        cdef double detF, detF0, detF_dt
        cdef np.ndarray[np.int64_t,ndim=1] udof, rdof
        cdef np.ndarray[double,ndim=1] ve, veE, ve0,ve0E, r0, r0E,N, wp, vec,drho_dt,gradRho, qhat
        cdef np.ndarray[double,ndim=2] gp, dNx, J, P
        cdef dict deform, deform0, M, Diff, RhoS

        udof = np.array([0, 1, 3, 4, 6, 7, 9, 10]);
        rdof = np.array([2, 5, 8, 11]);

        ve  = self.Dof(dc);
        veE = self.T.dot(ve);        # ve in element coordinates
        ve0  = self.Dof0(dc);
        ve0E = self.T.dot(ve0);        # ve0 in element coordinates
        r0  = self.Loc();
        r0E = self.T[np.ix_(udof,udof)].dot(r0);        # r0 in element coordinates


        drho_dt = ve - ve0;
        drho_dt = drho_dt[rdof]/dc.dt;

        gp = 0.577350269189626* \
             np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])

        wp = np.array([1.,1.,1.,1.]);
        vec = np.zeros(len(ve));
        ndim = 2;
        

        # Loop through integration points
        for ip in range(4):
            # r,s coordinates at this integration point
            r_ip = gp[ip,0];
            s_ip = gp[ip,1];

            # Compute deformation gradient
            deform = eh.LocalDeformation_Q2(r_ip,s_ip,r0E,veE[udof])
            N = deform['N']
            dNx = deform['dNx'].T
            J = deform['J'];
            detJ = np.linalg.det(J)
            Finv = np.linalg.inv(deform['F']);

            # Compute deformation gradient
            deform0 = eh.LocalDeformation_Q2(r_ip,s_ip,r0E,ve0E[udof])
            detF0 = np.linalg.det(deform0['F']);
            detF  = np.linalg.det(deform['F']);
            detF_dt = (detF-detF0)/dc.dt

            # Compute gradient of rho wrt reference X
            Rho = 0;
            dRhodt = 0;
            gradRho = np.zeros(ndim);
            for j in range(4):
                Rho    += N[j]*veE[rdof[j]];
                dRhodt += N[j]*drho_dt[j];
                for s in range(ndim):
                     gradRho[s] += dNx[s,j]*veE[rdof[j]];
            # Compute gradient of rho wrt current x
            gradRho = Finv.T.dot(gradRho);

            # Compute tangent matrix and updated deformation gradient
            M = self.property.Piola1Stiffness(deform['F'],eh.getFg(self,ip,dc,ndim), Rho, eh.getNc(self,ip,dc),0.,ndim);
            Diff = self.property.Diffusivity(deform['F'],Rho,gradRho,ndim);
            RhoS = self.property.SourceRho(deform['F'],Rho,gradRho,dc,ndim);
            qhat = Finv.dot(Diff['q'])
            P = M['P'];
            
            for i in range(4):
                eni = i*(ndim+1);        # +1 due to additional field
                
                # Displacement contribution
                for p in range(ndim):
                    for q in range(ndim):
                        vec[eni+p] += P[p,q]*dNx[q,i]*detJ*wp[ip];

                # Density contribution
                vec[eni+ndim]+= N[i]*(dRhodt+Rho*detF_dt/detF)*detJ*wp[ip];
                vec[eni+ndim]-= N[i]*RhoS['fRho']*detJ*wp[ip];
                for p in range(ndim):
                    vec[eni+ndim] += dNx[p,i]*qhat[p]*detJ*wp[ip]


        # Transform back into global x,y,z coordinates
        vec = self.T.T.dot(vec)

        return vec

    cpdef getVolume(self,DataContainer dc):
        cdef int ndim, ip,
        cdef double r_ip,s_ip, volume
        cdef np.ndarray[np.int64_t,ndim=1] udof
        cdef np.ndarray[double,ndim=1] wp, r
        cdef np.ndarray[double,ndim=2] gp
        cdef dict jac

        udof = np.array([0, 1, 3, 4, 6, 7, 9, 10]);
        r  = self.Loc()+self.Dof(dc)[udof];
        

        gp = 0.577350269189626* \
             np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])
        wp = np.array([1.,1.,1.,1.]);
        ndim = 2;

        volume = 0;

        # Loop through integration points
        for ip in range(4):
            # r,s coordinates at this integration point
            r_ip = gp[ip,0];
            s_ip = gp[ip,1];

            # Compute jacobian and update volume
            jac = eh.Jacobian_Q2(r_ip,s_ip,r)
            volume += np.linalg.det(jac['J'])

        return volume;

    cpdef getDensity(self,DataContainer dc):
        cdef np.ndarray[np.int64_t,ndim=1] rdof
        rdof = np.array([2, 5, 8, 11]);
        return np.mean(self.Dof(dc)[rdof]);