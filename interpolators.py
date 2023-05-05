import numpy as np
from numpy import linalg

class Cubic:
    
    # constructor function    
    def __init__(self, Xs = [], fs = [], Integrate = True, Extend = "Linear", Xmin = None, Xmax = None):
        self.Xs = Xs
        self.fs = fs
        self.dimension = len(Xs)
        
        # 0. Computes H
        self.H = [self.Xs[i]-self.Xs[i-1] for i in range(1,self.dimension)]
        
        # 1. Populates TH
        TH = [[0 for _ in range(self.dimension)] for _ in range(self.dimension)]
        for i in range(1,self.dimension-1):
            TH[i][i-1],TH[i][i],TH[i][i+1] = [self.H[i-1],2*(self.H[i-1]+self.H[i]),self.H[i]]
        # natural end conditions
        TH[0][0],TH[-1][-1] = [1,1]
        self.TH = TH
        
        # 2. Populates COB
        COB = [[0 for _ in range(self.dimension)] for _ in range(self.dimension)]
        for i in range(1,self.dimension-1):
            COB[i][i-1],COB[i][i],COB[i][i+1] = [3/self.H[i-1],-3/self.H[i-1]-3/self.H[i],3/self.H[i]]
        self.COB = np.matrix(COB)
        
        # 3. Solves C in basis of F
        CinF = np.matmul(linalg.inv(np.matrix(self.TH)),self.COB)
        self.CinF = CinF
        
        # 4. Solves B by b1 and propagation
        BinF = [[0 for _ in range(self.dimension)] for _ in range(self.dimension)]
        # b0 by Eq(18.21)
        BinF[0][0],BinF[0][1] = [-1/self.H[0],1/self.H[0]]
        BinF[0] += (-self.H[0]/3)*self.CinF[1]
        BinF[-1] = np.matrix(BinF[-1])
        # b next by Eq(18.23)
        i=1
        while i<self.dimension-1:
            BinF[i] = BinF[i-1] + self.H[i-1]*(self.CinF[i-1]+self.CinF[i])
            i += 1
        self.BinF = BinF
        
        # 5. Solves D by Eq(18.18)
        DinF = [[0 for _ in range(self.dimension)] for _ in range(self.dimension)]
        DinF[-1] = np.matrix(DinF[-1])
        i=0
        while i<self.dimension-1:
            DinF[i] = (self.CinF[i+1]-self.CinF[i])/(3*self.H[i])
            i += 1
        self.DinF = DinF
        
        # 6. Populates A
        AinF = np.zeros((self.dimension, self.dimension))
        np.fill_diagonal(AinF, 1)
        self.AinF = AinF
                
        # 7. Weight matrix of F for an integral of the interpolated segments
        W_inter = []
        for i in range(self.dimension-1):
            Wi = np.matrix(self.H[i]*self.AinF[i]\
                           +np.power(self.H[i],2)*self.BinF[i]/2\
                           +np.power(self.H[i],3)*self.CinF[i]/3\
                           +np.power(self.H[i],4)*self.DinF[i]/4)
            W_inter.append(np.matrix(Wi))
        
        # 8. Extrapolates and Updates Weights
        self.W_extra = W_inter
        if Extend in ["Linear","Natural"]:
            
            hmin = Xs[0]-Xmin
            
            Amin = hmin*AinF[0]
            self.AinF[0] += Amin
            Bmin = -np.power(hmin,2)*BinF[0]/2
            self.BinF[0] += Bmin
            
            self.Wmin = np.matrix(Amin + Bmin)
            
            hmax = Xmax-Xs[-1]
            
            Amax = hmax*AinF[-2]
            self.AinF[-2] += Amax
            Bmax = (np.power(hmax+self.H[-1],2)-np.power(self.H[-1],2))*BinF[-2]/2
            self.BinF[-2] += Bmax
            
            self.Wmax = np.matrix(Amax + Bmax)
            
            if Extend == "Natural":
                Cmin = 0
                Dmin = -np.power(hmin,4)*DinF[0]/4
                self.DinF[0] += Dmin
                self.Wmin += Dmin
                self.Wmin = np.matrix(self.Wmin)
                
                Cmax = (np.power(hmax+self.H[-1],3)-np.power(self.H[-1],3))*CinF[-2]/3
                self.CinF[-2] += Cmax
                self.Wmax += Cmax
                Dmax = (np.power(hmax+self.H[-1],4)-np.power(self.H[-1],4))*DinF[-2]/4
                self.DinF[-2] += Dmax
                self.Wmax += Dmax
                self.Wmax = np.matrix(self.Wmax)
            
            W_extra = [self.Wmin]
            for i in range(0,self.dimension-1):
                W_extra.append(self.W_extra[i])
            W_extra.append(self.Wmax)
            self.W_extra = W_extra
                
        # 9. Integrates in segments
        if Integrate == True:
            self.I_seg = np.matmul(self.W_extra,np.array(fs).T)
            self.I = sum(self.I_seg)
            self.Iint = sum(np.matmul(W_inter,np.array(fs).T))
            if Extend in ["Linear","Natural"]:
                self.Imin = sum(np.matmul(self.Wmin,np.array(fs).T))
                self.Imax = sum(np.matmul(self.Wmax,np.array(fs).T))
                
class Linear:

    # constructor function
    def __init__(self, Xs = [], fs = [], Integrate = True, Extend = "Linear", Xmin = None, Xmax = None):
        self.Xs = Xs
        self.fs = fs
        self.dimension = len(Xs)
        
        # 0. Computes H
        self.H = [self.Xs[i]-self.Xs[i-1] for i in range(1,self.dimension)]

        # 1. Populates A
        AinF = np.zeros((self.dimension, self.dimension))
        np.fill_diagonal(AinF, 1)
        self.AinF = AinF

        # 2. Solves B
        BinF = [[0 for _ in range(self.dimension)] for _ in range(self.dimension)]
        for i in range(1,self.dimension):
            BinF[i-1] += (self.AinF[i]-self.AinF[i-1])/self.H[i-1]
        BinF[-1] = np.matrix(BinF[-1])
        self.BinF = BinF

        # 3. Weight matrix of F for an integral of the interpolated segments
        W_inter = []
        for i in range(self.dimension-1):
            Wi = np.matrix(self.H[i]*self.AinF[i]\
                           +np.power(self.H[i],2)*self.BinF[i]/2)
            W_inter.append(np.matrix(Wi))

        # 4. Extrapolates and Updates Weights
        self.W_extra = W_inter
        if Extend in ["Linear"]:
            
            hmin = Xs[0]-Xmin
            
            Amin = hmin*AinF[0]
            self.AinF[0] += Amin
            Bmin = -np.power(hmin,2)*BinF[0]/2
            self.BinF[0] += Bmin
            
            self.Wmin = np.matrix(Amin + Bmin)
            
            hmax = Xmax-Xs[-1]
            
            Amax = hmax*AinF[-2]
            self.AinF[-2] += Amax
            Bmax = (np.power(hmax+self.H[-1],2)-np.power(self.H[-1],2))*BinF[-2]/2
            self.BinF[-2] += Bmax
            
            self.Wmax = np.matrix(Amax + Bmax)

            W_extra = [self.Wmin]
            for i in range(0,self.dimension-1):
                W_extra.append(self.W_extra[i])
            W_extra.append(self.Wmax)
            self.W_extra = W_extra
        
        # 5. Integrates in segments
        if Integrate == True:
            self.I_seg = np.matmul(self.W_extra,np.array(fs).T)
            self.I = sum(self.I_seg)
            self.Iint = sum(np.matmul(W_inter,np.array(fs).T))
            if Extend in ["Linear"]:
                self.Imin = sum(np.matmul(self.Wmin,np.array(fs).T))
                self.Imax = sum(np.matmul(self.Wmax,np.array(fs).T))
