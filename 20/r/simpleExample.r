source("https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/r/y20_subroutine.r")
#source("y20_subroutine.r")  #In case you have it local on your computer

#RhpcBLASctl::blas_set_num_threads(1) # Single core

theta <- c(0.51,5.19,0.13,0.1,0.5,0.,1.,1.,0.99,0.,0.,0.,0.,0.,0.163,0.,-0.,0.,0.,0.,0.,0.158,-0.002,0.17,-0.005,0.067,-0.,-1.44,-2.0,-6.9,0.0042,0.0015,-2.55,1.24,0.25)
init <- rep(0,5)
infall <- c(0.5,0.1,0.1,0.2,0) #0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
time <- 1   #Time to run
avgT <- c(-2., 2., 5., 10., 15., 20., 22., 21., 16., 11., 6., 3.)  #Temp month average [C]
sumP <- 600 #Precip annual summ [mm]
diam <- 2   #size [cm]
leach <- 0  #Leaching

yasso20$setTheta(theta)
yasso20$setA(avgT, sumP, diam, leach)

for(year in 0:9) {
  res <- yasso20$getNextTimestep(init, infall, time)
  print(c(year,res));
  init <- res
}
res <- yasso20$getSpin(infall)
print(res)

system.time({
  n <- 1000000
  result <- rep(0, 5)
  yasso20$setTheta(theta)
  for(i in 1:n) {
    yasso20$setA(avgT + i / n, sumP, diam, leach)
    result  <- yasso20$getNextTimestep(result, infall * i / n, time)
  }
})
# Single Core
#       User      System verstrichen 
#    113.378       0.023     113.483 

# Multiple Core
#      User      System verstrichen 
#    220.483     454.083     134.808 

result
#[1] 1.4053736 0.1492327 0.2233364 3.0202410 6.6699129
