## ---- packages ----
source(here::here('00_packages/packages.R'))

## ---- simulation ----
# On simulation of tempered stable random variates


## --- Algorithm 0 ----

# simulation parameters
b <- 0.1
alpha <- 1.2 # alpha betwenn 0 and 1
a <- 1 # delta_+ LV 215
# Levy matters IV -> page 215
# beta -> 1; Subordinator; absolut skewness; 
# delta -> 0; location
# gamma -> sigma^beta in LV and there delta_+ equals -> a in the simulation paper

sigma_alpha <- (1/alpha) * gamma(1- alpha) * a * cos((alpha*pi)/2)
gamma <- sigma_alpha^(1/alpha)

Y <- NULL
N <- 100
i <-0 

while (length(Y) < N) {
  U <- runif(1)
  V <- stabledist::rstable(1, alpha = alpha, beta = 1, gamma = gamma, delta = 0)
  if(U <= exp(-b*V)) Y <- c(Y, V)
  i <- i + 1
  print(paste0('Y = ', length(Y), '; ', 'tries = ', i))
}


## --- Algorithm 2 ----

## simulation parameters
c <- 0.06 # > 0
alpha <- 1.2 # alpha betwenn 0 and 1
b <- 1
a <- 1 # delta_+ LV 215
delta_a <- 0.001
# Levy matters IV -> page 215
# beta -> 1; Subordinator; absolut skewness; 
# delta -> 0; location
# gamma -> sigma^beta in LV and there delta_+ equals -> a in the simulation paper

sigma_alpha <- (1/alpha) * gamma(1- alpha) * (alpha*delta_a) * cos((alpha*pi)/2)
gamma <- sigma_alpha^(1/alpha)

Y_1 <- NULL
Y_2 <- NULL
N <- 100
i <-0 

# berechnungen 

U <- runif(1)
V <- stabledist::rstable(1, alpha = alpha, beta = 1, gamma = gamma, delta = 0)
if(U <= exp(-b*(V + c))) Y_2 <- V - delta_a * gamma(1 - alpha) *a * b^(alpha *-1)




## loop 
while (length(Y) < N) {
  i <- i + 1
  print(paste0('Y = ', length(Y), '; ', 'tries = ', i))
}









