# -*- coding: utf-8 -*-\
# Generalized maxwell model, with sigma_v = 2*eta*d_v
from Property cimport *
from math import *
import numpy as np
cimport numpy as np

cdef class PMAXWELL(Property):
    def __init__(self,material0,material1,propProp):
        super(PMAXWELL,self).__init__('PMAXWELL');
        
        self.eta     = propProp['eta'];
        self.material0 = material0;
        self.material1 = material1;

    def __str__(self):
        return super(PMAXWELL,self).__str__() + "\t eta = " + str(self.eta) + "\t material0 = " + str(self.material0.localID) + "\t material1 = " + str(self.material1.localID)

    cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, double dt, int ndim):
        cdef int i,j,k,l,p,q,r,s,t
        cdef double J
        cdef np.ndarray[double, ndim=2] Av,Bv,Avinv, Fe, Fv,Finv,Feinv,Fvinv,P, P0,Pe,P1
        cdef np.ndarray[double, ndim=4] A, A0, A1, dPdFv, dFvdF, dBdFe
        cdef dict M0, M1

        if(dt>0):
            # Update fv
            #self.UpdateFg(F,Fg0,dt,ndim)
            #Fv = Fg0

            Fv, Av, Bv = self.UpdateFg(F,Fg0,dt,ndim)
        else:
            Fv = Fg0
        Fe = np.linalg.solve(Fv.T,F.T)
        Feinv = np.linalg.inv(Fe)
        Fvinv = np.linalg.inv(Fv)

        # Contribution from different materials
        M0 = self.material0.Piola1Stiffness(F,ndim);
        M1 = self.material1.Piola1Stiffness(Fe,ndim);
        
        # Compute stresses
        P0 = M0['Pe']
        Pe = M1['Pe']
        P1 = M1['Pe'].dot(Fvinv.T)
        P = P0+P1


        if(dt>0):
            A0 = M0['Ae'];
            A1 = M1['Ae'];

            Finv  = np.linalg.inv(F)

            # Compute dP_ij/dFv_rs
            dPdFv = np.zeros((ndim,ndim,ndim,ndim))
            for i in range(ndim):
                for j in range(ndim):
                    for r in range(ndim):
                        for s in range(ndim):
                            dPdFv[i,j,r,s] = -P1[i,s]*Fvinv[j,r]
                            for p in range(ndim):
                                for q in range(ndim):
                                    for t in range(ndim):
                                        dPdFv[i,j,r,s] -=Fvinv[j,p]*A1[i,p,q,t]*Fvinv[s,t]*Fe[q,r]

            """
            # Compute dB_ij/dFe_kl
            dBdFe = np.zeros((ndim,ndim,ndim,ndim))
            for i in range(ndim):
                for j in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim):
                            for p in range(ndim):
                                dBdFe[i,j,k,l] += Feinv[i,p]*Pe[p,l]*Fe[k,j]
                                for q in range(ndim):
                                    dBdFe[i,j,k,l] += Feinv[i,p]*Pe[p,q]*Fe[k,q]*(j==l)
                                    for r in range(ndim):
                                        dBdFe[i,j,k,l] += -Feinv[p,l]*Feinv[k,l]*Pe[p,q]*Fe[r,q]*Fe[r,j] + \
                                                           Feinv[i,p]*A0[p,q,k,l]*Fe[r,q]*Fe[r,j]
                            dBdFe[i,j,k,l]/=(2.*self.eta)
            """
            # Compute dFv_rs/dF_kl
            J = np.linalg.det(F)
            Avinv = np.linalg.inv(Av)
            dFvdF = np.zeros((ndim,ndim,ndim,ndim))
            for r in range(ndim):
                for s in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim):
                            for i in range(ndim):
                                for p in range(ndim):
                                    dFvdF[r,s,k,l] +=dt/J*Avinv[r,i]*Bv[i,p]*\
                                           (dFvdF[p,s,k,l]-Fv[p,s]*Finv[l,k])
                                    """
                                    for q in range(ndim):
                                        dFvdF[r,s,k,l] +=dt/J*Avinv[r,i]*dBdFe[i,p,r,q]*\
                                           (r==k-Fe[r,k])*Fvinv[l,q]*Fv[p,s]
                                    """


            # Compute A_ijkl
            A = A0;
            for i in range(ndim):
                for j in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim):
                            for p in range(ndim):
                                for q in range(ndim):
                                    A[i,j,k,l] += Fvinv[j,p]*A1[i,p,k,q]*Fvinv[l,q] +\
                                                  dPdFv[i,j,p,q]*dFvdF[p,q,k,l]

            
            return {'P':P,'A':A,'Fg':Fv}
        else:
            return{'P':P}


    cpdef UpdateFg(self, np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fv0,double dt,int ndim):
        cdef np.ndarray[double,ndim=2] Fe, Fv, A, B, R
        cdef np.ndarray[double,ndim=2] Fe1,Fv1,A1,B1,R1,dFv1,delta
        cdef np.ndarray[double, ndim=4] Ae, H
        cdef dict Me

        cdef double dR, J, normF,normF0 
        cdef int i,j,k,l,it,itt

        Fv = Fv0
        J = np.linalg.det(F)

        delta = np.eye(ndim)
        dR = 1.e-2

        
        it = 0;
        epsilon = 1;
        while(epsilon>1e-12 and it<50):
            it+=1;

            Fe = np.linalg.solve(Fv.T,F.T)
            Me = self.material1.Piola1Stiffness(Fe,ndim);
            B  = np.linalg.inv(Fe).dot(Me['Pe']).dot(Fe.T).dot(Fe)/(2.*self.eta)
            A  = delta-B*(dt/J)

            R = A.dot(Fv)-Fv0

            epsilon =  np.linalg.norm(R)


            # Compute H_ijkl = dFv_ij/dR_kl
            H      = np.zeros((ndim,ndim,ndim,ndim))
            Hcheck = np.zeros((ndim,ndim,ndim,ndim))
            for k in range(ndim):
                for l in range(ndim):
                    A1 = A 
                    R1 = R+dR*np.outer(delta[:,k],delta[:,l])

                    Fv1  = np.linalg.solve(A1,R1+Fv)

                    normF = np.inf 
                    norm  = 0;
                    itt = 0;
                    while (np.abs(normF-norm0)>1e-12 and itt<50):
                        itt+=1 
                        norm0 = normF

                        Fe1 = np.linalg.solve(Fv1.T,F.T)
                        Me = self.material1.Piola1Stiffness(Fe1,ndim);
                        B1  = np.linalg.inv(Fe1).dot(Me['Pe']).dot(Fe1.T).dot(Fe1)/(2.*self.eta)
                        A1  = delta-B1*(dt/J)
                        Fv1  = np.linalg.solve(A1,R1+Fv)
                        dFv1 = (Fv1-Fv)/dR

                        normF = np.linalg.norm(dFv1)
                        """
                        print np.abs(normF-norm0)
                    print '\n'
                    """
            
                    H[:,:,k,l] = dFv1


            Ainv = np.linalg.inv(A)
            for i in range(ndim):
                for j in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim):
                            Hcheck[i,j,k,l] = Ainv[i,k]*(j==l)
            """
            print H
            print '\n'
            print Hcheck
            """

            dFv = np.zeros((ndim,ndim))
            for i in range(ndim):
                for j in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim): 
                            dFv[i,j] -= H[i,j,k,l]*R[k,l]

            
            # Update Fv
            #Fv  = Fv-np.linalg.solve(A,R)
            Fv  = Fv+dFv

            """
            print epsilon
        print '\n'
        """

        

        return (Fv,A,B)
