## ---- packages ----
source(here::here('00_packages/packages.R'))

## ---- Gamma Process simulation ----
# From Stoch Sim 
#####Exercise 11.2#####

#there are two possibilities to construct the variance gamma process

#(1) subordination of a Brownian motion
lambda <- 10
alpha <- 2
beta <- 1

tn <- 1/1000

Gincr <- rgamma(1000,shape = lambda*tn, rate = (alpha^2-beta^2)/2)
Tt <- cumsum(Gincr)
plot(seq(0,1,tn),c(0,Tt),main = "Gamma subordinator",xlab = "t",ylab = "T(t)")

VGincr <- rnorm(1000, beta*Gincr, sqrt(Gincr))
Xt <- cumsum(VGincr)
plot(seq(0,1,tn),c(0,Xt),main = "Variance Gamma process",xlab = "t",ylab = "X(t)")


#(2) difference of two gamma processes
C <- lambda
G <- alpha+beta
M <- alpha-beta

Gincrp <- rgamma(1000,shape = C*tn, rate = M)
Gincrm <- rgamma(1000,shape = C*tn, rate = G)

Xp <- cumsum(Gincrp)
Xm <- cumsum(Gincrm)

Xt <- Xp-Xm
plot(seq(0,1,tn),c(0,Xt),main = "variance Gamma process",xlab = "t",ylab = "X(t)")

















