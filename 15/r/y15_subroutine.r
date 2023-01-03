library(Matrix)

yasso15 <- new.env()

yasso15$setTheta <- function(theta_) {
  yasso15$theta <- theta_
  tabs <- c(1:4, 32, 35)
  yasso15$theta[tabs] = -abs(yasso15$theta[tabs])
}

yasso15$setA <- function(avgT, sumP, ampT, diam, leach) {
  yasso15$A <- matrix(0, 5, 5)

  m3 <- sumP/1000
#temperature annual cycle approximation
  m0 <- (1./sqrt(2.)-1.)/pi
  m1 <- 1./sqrt(2.)/pi
  m2 <- (1.-1./sqrt(2.))/pi
  clim <- 4. * ampT
  te <- 1:4
  te[1] = avgT + clim * m0
  te[2] = avgT - clim * m1
  te[3] = avgT + clim * m2
  te[4] = avgT + clim * m1

#Average temperature dependence
  tem = sum(exp(yasso15$theta[22]*te+yasso15$theta[23]*te^2)) / 4
  temN = sum(exp(yasso15$theta[24]*te+yasso15$theta[25]*te^2)) / 4
  temH = sum(exp(yasso15$theta[26]*te+yasso15$theta[27]*te^2)) / 4

#Precipitation dependence
  tem = tem * (1.-exp(yasso15$theta[28] * m3));
  temN = temN * (1.-exp(yasso15$theta[29] * m3));
  temH = temH * (1.-exp(yasso15$theta[30] * m3));

#Size class dependence -- no effect if d == 0.0
  size_dep <- 1
  if(diam > 0.) {size_dep <- min(1., (1. + yasso15$theta[33]*diam + yasso15$theta[34]*diam^2)^yasso15$theta[35])}
#Calculating matrix a (will work ok despite the sign of alphas)
  yasso15$A[1+(0:2*6)] = yasso15$theta[1:3]*tem*size_dep
  yasso15$A[4,4] = yasso15$theta[4]*temN*size_dep
  yasso15$A[5,5] = yasso15$theta[32]*temH #no size effect in humus
  dAbs <- abs(yasso15$A[1+(0:3*6)])
  idx <- 5
  for(i in 0:3) {
    for(j in 0:3) {
      if(i!=j) {
        yasso15$A[1+j*5+i] = yasso15$theta[idx] * dAbs[1+j];
        idx <- idx+1
      }
    }
  }
#mass flows AWEN -> H (size effect is present here)
  yasso15$A[1:4*5] = yasso15$theta[31] * dAbs
#Leaching (no leaching for humus) 
  if(leach < 0.) {yasso15$A[1+6*0:3] = yasso15$A[1+6*0:3] + leach * m3}
}

yasso15$getSpin <- function(infall) {
  solve(yasso15$A, infall) * -1
}

yasso15$getNextTimestep <- function(init, infall, time) {
  z1 <- yasso15$A %*% init
  z1 <- z1 + infall
  At <- yasso15$A * time
  z2 <- Matrix::expm(At) %*% z1
#Alternatives 
#library(expm)
#z2 <- expAtv(At,z1)$eAtv
#
#library(expoRkit)
#z2 <- expv(At, z1)
  z2 <- z2 - infall
  as.vector(solve(yasso15$A, z2))
}
