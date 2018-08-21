# -*- coding: utf-8 -*-\
# Anisotropic prescribed growth: Fg = I + (theta-1)*e1*e1
from Property cimport *
from math import *
from copy import deepcopy
import numpy as np
cimport numpy as np

cdef class PSOLID71(Property):
    def __init__(self,material,propProp):
        super(PSOLID71,self).__init__('PSOLID71');
        
        self.Grho               = propProp['Grho'];
        self.advSpeed           = propProp['advSpeed']
        self.kth1               = propProp['kth1'];
        self.kth2               = propProp['kth2'];
        self.alpha              = propProp['alpha'];
        self.amplRhoSource      = propProp['amplRhoSource'];
        try:
            self.amplAdvThreshold = propProp['amplAdvThreshold'];
        except:
            self.amplAdvThreshold = None
        self.material           = material;

    def __str__(self):
        return super(PSOLID71,self).__str__() + "\t Grho = " + str(self.Grho)  \
                                              + "\t advSpeed = " + str(self.advSpeed)  \
                                              + "\t alpha = " + str(self.alpha)  \
                                              + "\t kth1 = " + str(self.kth1)  \
                                              + "\t kth2 = " + str(self.kth2) \
                                              + "\t material = " + str(self.material.localID)

    cpdef Piola1Stiffness(self,np.ndarray[double,ndim=2] F, np.ndarray[double,ndim=2] Fg0, \
                          double Rho, double nc0, double dt, int ndim):
        cdef int i,j,k,l
        cdef double th1, th2, thRat, th1Coeff, th2Coeff, Jg
        cdef np.ndarray[double, ndim=1] e1
        cdef np.ndarray[double, ndim=2] Fe, Fg, Fginv, P, dPdRho
        cdef np.ndarray[double, ndim=4] A, Ae

        # Update nc, th1, th2
        if(dt>0):
            nc = self.UpdateNc(nc0,dt)
        else:
            nc = nc0
        th1, dth1dRho = self.UpdateTh1(nc,Rho);
        th2, dth2dRho = self.UpdateTh2(nc,Rho);


        if ndim==3:
            Fg    = np.diag(np.array([th2,th1,th1]));
            Fginv = np.diag(1./np.array([th2,th1,th1]));
            Jg    = th2*th1*th1
        else:
            Fg    = np.diag(np.array([th2,th1]));
            Fginv = np.diag(1./np.array([th2,th1]));
            Jg    = th2*th1
        #Jg=1
        Fe = F.dot(Fginv);
        
        e1 = np.zeros(ndim); e1[0] = 1.;
        thRat = (th1-th2)/th2;

        # Compute elastic stress and elastic tangent
        Me = self.material.Piola1Stiffness(Fe,ndim);

        # Compute P_ij = (Pe_ij + thRat*Pe_i1*e1_j)/th1
        P = Jg*Me['Pe']
        for i in range(ndim):
            P[i,0]+= thRat*P[i,0];
        P/=th1;

        if dt>0:
            # A_ijkl = dP_ij/dF_ij
            Ae = Jg*Me['Ae']
            A = np.zeros((ndim,ndim,ndim,ndim))
            for i in range(ndim):
                for j in range(ndim):
                    for k in range(ndim):
                        for l in range(ndim):
                            A[i,j,k,l] = Ae[i,j,k,l] + thRat* (Ae[i,j,k,0]*e1[l] + Ae[i,0,k,l]*e1[j]) \
                                                     + thRat**2 * Ae[i,0,k,0]*e1[j]*e1[l]
            A/=(th1**2);

            # dPdRho_ij = dP_ij/dRho
            Fgid = np.diag(Fginv);
            if ndim==3:
                dFgd = np.array([dth2dRho,dth1dRho,dth1dRho]);
            else:
                dFgd = np.array([dth2dRho,dth1dRho]);
              
            dPdRho = np.zeros((ndim,ndim))
            for i in range(ndim):
                for j in range(ndim):
                    dPdRho[i,j] -= P[i,j]*Fgid[j]*dFgd[j]
                    dPdRho[i,j] += P[i,j]*np.dot(Fgid,dFgd)
                    for k in range(ndim):
                        for l in range(ndim):
                            dPdRho[i,j] -=  Fgid[j]*Fgid[l]*dFgd[l]*Ae[i,j,k,l]*Fe[k,l]
            
           

            """
            th1Coeff = -dth1dRho/th1/th1
            th2Coeff = -dth2dRho/th2/th2
            dPdRho = deepcopy(Me['Pe'])
            dPdRho[0,0]*=th2Coeff
            dPdRho[1,1]*=th1Coeff
            if ndim==3:
                dPdRho[2,2]*=th1Coeff
            """
        else:
            A = np.zeros((ndim,ndim,ndim,ndim))
            dPdRho = np.zeros((ndim,ndim))

        return {'P':P,'A':A,'dPdRho':dPdRho,'Fg':Fg,'nc':nc}


    cpdef Diffusivity(self,np.ndarray[double,ndim=2] F, \
                          double Rho, np.ndarray[double,ndim=1] gradRho, int ndim):

        cdef double stretch_v
        cdef np.ndarray[double, ndim=1] q, D1, v,v_tilde, advFlux, advTang
        cdef np.ndarray[double, ndim=2] D2
        cdef np.ndarray[double, ndim=3] advTangF, D3


        # Compute advection speed velocity and tangents
        v_tilde = F[:,0]
        stretch_v = np.linalg.norm(v_tilde)
        v = self.advSpeed*v_tilde/stretch_v

        advFlux = -v*Rho
        advTang = -v
        advTangF = np.zeros((ndim,ndim,ndim))
        for i in range(ndim):
            for k in range(ndim):
                advTangF[i,k,0] = -Rho*self.advSpeed* \
                                   ( (i==k)/stretch_v \
                                     - v_tilde[i]*v_tilde[k]/(stretch_v**3) )

        if not self.amplAdvThreshold == None:
            advFlux *= self.amplAdvThreshold.Get(Rho)
            advTang *= (self.amplAdvThreshold.Get(Rho) + \
                                  Rho*self.amplAdvThreshold.GetDerivative(Rho))
            advTangF*= self.amplAdvThreshold.Get(Rho)



        """
        v = v_tilde*self.advSpeed
        if self.amplAdvThreshold == None:
            advFlux = -v*Rho
            advTang = -v
            advTangF = np.zeros((ndim,ndim,ndim));
        else:
            advFlux = -v*Rho*self.amplAdvThreshold.Get(Rho)
            advTang = -v*(self.amplAdvThreshold.Get(Rho) + \
                                  Rho*self.amplAdvThreshold.GetDerivative(Rho))
            advTangF = np.zeros((ndim,ndim,ndim));
        """

        q = self.material.D*gradRho + advFlux;
        D1 = advTang;                           # = dq/dRho       
        D2 = self.material.D*np.eye(ndim);      # = dq/dGradRho
        D3 = advTangF;        # = dq/dF

        """
        print "self.advSpeed: ", self.advSpeed
        print "v_tilde:       ", v_tilde
        print "stretch_v:     ", stretch_v
        print "v:             ", v
        print "norm(v):       ", np.linalg.norm(v)
        print "\n\n\n"
        """
        
        return {'q':q,'D1':D1,'D2':D2,'D3':D3};

    cpdef SourceRho(self,np.ndarray[double,ndim=2] F, \
                    double Rho, np.ndarray[double,ndim=1] gradRho, \
                    DataContainer dc, int ndim):

        cdef double c1, beta, fRho, dfdRho,J;
        cdef np.ndarray[double, ndim=2] dfdF, Finv;

        """
        c1 = self.Grho*self.amplRhoSource.Get(dc.time);
        beta = .5;
        fRho = c1*(beta-Rho);
        dfdRho = -c1;
        dfdF = np.zeros((ndim,ndim));
        """
        c1 = self.Grho*self.amplRhoSource.Get(dc.time);
        fRho = c1;
        dfdRho = 0;
        dfdF = np.zeros((ndim,ndim))

        return {'fRho':fRho,'dfdRho':dfdRho,'dfdF':dfdF};

    cpdef UpdateNc(self,double nc0, double dt):
        return nc0+0*dt;

    cpdef UpdateTh1(self,double nc, double Rho):
        return (            (1.+self.kth1*Rho)**self.alpha, 
                 self.alpha*(1.+self.kth1*Rho)**(self.alpha-1.)*self.kth1);

    cpdef UpdateTh2(self,double nc, double Rho):
        return (            (1.+self.kth2*Rho)**self.alpha, 
                 self.alpha*(1.+self.kth2*Rho)**(self.alpha-1.)*self.kth2);
