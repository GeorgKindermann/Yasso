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
    infall2 = infall * i / n
    avgT2 = avgT + i / n
    y20.setA(avgT2, sumP, diam, leach)
    result = y20.getNextTimestep(result, infall, time)

print(timeit.default_timer() - start_time)
#307.1552226519998  when using 5 threads
#108.64224406100038

result
#array([1.40538097, 0.14923343, 0.22333682, 3.02026497, 6.67175339])


