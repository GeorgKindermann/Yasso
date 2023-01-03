import numpy as np
from scipy.linalg import expm, solve

class yasso20:
    theta = np.array([0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25])
    A = np.zeros((5, 5))
    def setTheta(self, theta):
        self.theta = theta
        tabs = np.r_[0:4,[31,34]]
        self.theta[tabs] = -abs(self.theta[tabs])
    def setA(self, avgT, sumP, diam, leach):
        m3 = sumP/1000
        tem = np.sum(np.exp(self.theta[21]*avgT + self.theta[22]*avgT**2)) / 12
        temN = np.sum(np.exp(self.theta[23]*avgT + self.theta[24]*avgT**2)) / 12
        temH = np.sum(np.exp(self.theta[25]*avgT + self.theta[26]*avgT**2)) / 12
        #Precipitation dependence
        tem = tem * (1.-np.exp(self.theta[27] * m3));
        temN = temN * (1.-np.exp(self.theta[28] * m3));
        temH = temH * (1.-np.exp(self.theta[29] * m3));
        #Size class dependence -- no effect if d == 0.0
        size_dep = np.min([1., (1. + self.theta[32]*diam + self.theta[33]*diam**2)**self.theta[34]]) if diam > 0 else 1
        #Calculating matrix A (will work ok despite the sign of alphas)
        self.A[[0, 1, 2], [0, 1, 2]] = self.theta[0:3]*tem*size_dep
        self.A[3,3] = self.theta[3]*temN*size_dep
        self.A[4,4] = self.theta[31]*temH #no size effect in humus
        dAbs = np.abs(self.A[[0, 1, 2, 3], [0, 1, 2, 3]])
        idx = 4
        for i in range(4):
            for j in range(4):
                if i != j:
                    self.A[i,j] = self.theta[idx] * dAbs[j]
                    idx += 1
        #mass flows AWEN -> H (size effect is present here)
        self.A[4,:4] = self.theta[30] * dAbs
        #Leaching (no leaching for humus)
        if leach < 0.:
            self.A[[0, 1, 2, 3], [0, 1, 2, 3]] += leach * m3
    def getSpin(self,infall):
        return(-solve(self.A, infall))
    def getNextTimestep(self,init, infall, time):
        return(solve(self.A, expm(self.A * time).dot(self.A.dot(init) + infall)- infall))

