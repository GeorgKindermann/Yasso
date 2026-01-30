#To run it only on one core like all the other methods
import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import y20_subroutine as YASSO

import numpy as np

theta = np.array([0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25])
init = np.zeros(5)
infall = np.array([0.5,0.1,0.1,0.2,0])  #0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
time = 1   #Time to run
avgT = np.array([-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.])  #Temp month average [C]
sumP = 600 #Precip annual summ [mm]
diam = 2   #size [cm]
leach = 0  #Leaching


y20 = YASSO.yasso20()
y20.setTheta(theta)
y20.setA(avgT, sumP, diam, leach)

for year in range(10):
    res = y20.getNextTimestep(init, infall, time)
    print(year, res)
    init = res

y20.getSpin(infall)


import timeit
start_time = timeit.default_timer()
n = 1000000
infall2 = infall
result = np.zeros(5)
y20.setTheta(theta)
for i in range(n):
    y20.setA(avgT + i / n, sumP, diam, leach)
    result = y20.getNextTimestep(result, infall * i / n, time)

print(timeit.default_timer() - start_time)
# Single Core: 70.15
# Multicore: 87.00

result
# 1.40537233, 0.14923256, 0.2233362 , 3.02023804, 6.66990664]
