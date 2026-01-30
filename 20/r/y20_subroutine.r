library(Matrix)

yasso20 <- new.env()

yasso20$setTheta <- function(theta_) {
  yasso20$theta <- theta_
  tabs <- c(1:4, 32, 35)
  yasso20$theta[tabs] = -abs(yasso20$theta[tabs])
}

yasso20$Vidx <- c(6, 11, 16, 2, 12, 17, 3, 8, 18, 4, 9, 14)
yasso20$Aidx <- c(2, 3, 4, 1, 3, 4, 1, 2, 4, 1, 2, 3)

yasso20$setA <- function(avgT, sumP, diam, leach) {
  yasso20$A <- matrix(0, 5, 5)

  m3 <- sumP/1000
  tem = sum(exp(yasso20$theta[22]*avgT+yasso20$theta[23]*avgT^2)) / 12
  temN = sum(exp(yasso20$theta[24]*avgT+yasso20$theta[25]*avgT^2)) / 12
  temH = sum(exp(yasso20$theta[26]*avgT+yasso20$theta[27]*avgT^2)) / 12
  
#Precipitation dependence
  tem = tem * (1.-exp(yasso20$theta[28] * m3));
  temN = temN * (1.-exp(yasso20$theta[29] * m3));
  temH = temH * (1.-exp(yasso20$theta[30] * m3));

#Size class dependence -- no effect if d == 0.0
  size_dep <- 1
  if(diam > 0.) {size_dep <- min(1., (1. + yasso20$theta[33]*diam + yasso20$theta[34]*diam^2)^yasso20$theta[35])}
#Calculating matrix a (will work ok despite the sign of alphas)
  yasso20$A[1+(0:2*6)] = yasso20$theta[1:3]*tem*size_dep
  yasso20$A[4,4] = yasso20$theta[4]*temN*size_dep
  yasso20$A[5,5] = yasso20$theta[32]*temH #no size effect in humus
  dAbs <- abs(yasso20$A[1+(0:3*6)])
  #idx <- 5
  #for(i in 0:3) {
  #  for(j in 0:3) {
  #    if(i!=j) {
  #      yasso20$A[1+j*5+i] = yasso20$theta[idx] * dAbs[1+j];
  #      idx <- idx+1
  #    }
  #  }
  #}
  yasso20$A[yasso20$Vidx] <- yasso20$theta[5:16] * dAbs[yasso20$Aidx]
#mass flows AWEN -> H (size effect is present here)
  yasso20$A[1:4*5] = yasso20$theta[31] * dAbs
#Leaching (no leaching for humus) 
  if(leach < 0.) {yasso20$A[1+6*0:3] = yasso20$A[1+6*0:3] + leach * m3}
}

yasso20$getSpin <- function(infall) {
  solve(yasso20$A, infall) * -1
}

yasso20$getNextTimestep <- function(init, infall, time) {
  z1 <- yasso20$A %*% init
  z1 <- z1 + infall
  At <- yasso20$A * time
  z2 <- Matrix::expm(At) %*% z1
#Alternatives 
#library(expm)
#z2 <- expAtv(At,z1)$eAtv
#
#library(expoRkit)
#z2 <- expv(At, z1)
  z2 <- z2 - infall
  as.vector(solve(yasso20$A, z2))
}
