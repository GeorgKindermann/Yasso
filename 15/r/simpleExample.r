source("https://raw.githubusercontent.com/GeorgKindermann/Yasso/master/15/r/y15_subroutine.r")
#source("y15_subroutine.r")  #In case you have it local on your computer


theta <- c(0.49,4.9,0.24,0.095,0.44,0.25,0.92,0.99,0.084,0.011,0.00061,0.00048,0.066,0.00077,0.1,0.65,-0.15,-0.02,-0.92,-0.0004,-0.00017,0.091,-0.00021,0.049,-7.90E-05,0.035,-0.00021,-1.8,-1.2,-13,0.0046,0.0013,-0.44,1.3,0.26)
init <- rep(0,5)
infall <- c(0.5,0.1,0.1,0.2,0) #0..Acid soluble (Celluloses), 1..Water sol.(sugar), 2..Ethanol sol.(waxes), 3..Insoluble (lignin), 4..Humus
time <- 1   #Time to run
avgT <- 10  #Temp annual average [C]
sumP <- 600 #Precip annual summ [mm]
ampT <- 12  #Amplitude (max. difference of month averages / 2) [C]
diam <- 2   #size [cm]
leach <- 0  #Leaching

yasso15$setTheta(theta)
yasso15$setA(avgT, sumP, ampT, diam, leach)

for(year in 0:9) {
  res <- yasso15$getNextTimestep(init, infall, time)
  print(c(year,res));
  init <- res
}
res <- yasso15$getSpin(infall)
print(res)


system.time({
  n <- 1000000
  infall2 <- infall
  result <- rep(0, 5)
  yasso15$setTheta(theta)
  for(i in 1:n) {
    infall2  <- infall * i / n
    avgT2 <- avgT + i / n
    yasso15$setA(avgT2, sumP, ampT, diam, leach)
    result  <- yasso15$getNextTimestep(result, infall, time)
  }
})
#      User      System verstrichen 
#   1527.777       0.092    1528.432 
result
#[1]  2.6132458  0.2754559  0.3923254  8.3531430 10.5658222
