import numpy as np
from scipy.linalg import solve

class yasso20:
    def __init__(self, theta=None):
        # Initialisierung der Parameter
        if theta is None:
            self.theta = np.array([0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25])
        else:
            self.theta = theta.copy()
        self.setTheta(self.theta)
        self.A = np.zeros((5, 5))
    def setTheta(self, theta):
        self.theta = theta
        tabs = np.r_[0:4, [31, 34]]
        self.theta[tabs] = -np.abs(self.theta[tabs])
    def setA(self, avgT, sumP, diam, leach):
        th, m3 = self.theta, sumP / 1000.0
        tem = (np.exp(th[21]*avgT + th[22]*avgT**2).sum() / 12.0) * (1.0 - np.exp(th[27] * m3))
        temN = (np.exp(th[23]*avgT + th[24]*avgT**2).sum() / 12.0) * (1.0 - np.exp(th[28] * m3))
        temH = (np.exp(th[25]*avgT + th[26]*avgT**2).sum() / 12.0) * (1.0 - np.exp(th[29] * m3))
        size_dep = min(1.0, (1.0 + th[32]*diam + th[33]*diam**2)**th[34]) if diam > 0 else 1.0
        self.A.fill(0.0)
        d = np.zeros(4)
        d[0:3] = th[0:3] * (tem * size_dep)
        d[3] = th[3] * (temN * size_dep)
        if leach < 0: 
            d += leach * m3
        for i in range(4): 
            self.A[i, i] = d[i]
        self.A[4, 4] = th[31] * temH
        dAbs, idx = np.abs(d), 4
        for i in range(4):
            for j in range(4):
                if i != j:
                    self.A[i, j] = th[idx] * dAbs[j]
                    idx += 1
        self.A[4, :4] = th[30] * dAbs
    def getSpin(self, infall):
        return -solve(self.A, infall)
    def getNextTimestep(self, init, infall, time):
        vals, vecs = np.linalg.eig(self.A)
        vec_inv = np.linalg.solve(vecs, np.eye(5))
        expLt = np.exp(vals * time)
        intLt = np.where(np.abs(vals) > 1e-9, (expLt - 1.0) / vals, time)
        return (vecs @ (expLt * (vec_inv @ init) + intLt * (vec_inv @ infall))).real
