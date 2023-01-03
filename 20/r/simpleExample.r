source("https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/20/r/y20_subroutine.r")
#source("y20_subroutine.r")  #In case you have it local on your computer


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
  infall2 <- infall
  result <- rep(0, 5)
  yasso20$setTheta(theta)
  for(i in 1:n) {
    infall2  <- infall * i / n
    avgT2 <- avgT + i / n
    yasso20$setA(avgT2, sumP, diam, leach)
    result  <- yasso20$getNextTimestep(result, infall, time)
  }
})
#       User      System verstrichen 
#    224.090       0.000     224.125 

result
#[1] 1.4053808 0.1492334 0.2233368 3.0202649 6.6717529
