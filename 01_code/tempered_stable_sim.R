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
Delta <- 1
Delta_a <- a * Delta #0.001
# Levy matters IV -> page 215
# beta -> 1; Subordinator; absolut skewness; 
# delta -> 0; location
# gamma -> sigma^beta in LV and there delta_+ equals -> a in the simulation paper

sigma <- ( (1/alpha) * gamma(1- alpha) * (a) * cos((alpha*pi)/2) )^(1/alpha)

# gamma transformieren für groß-delta t^1/alpha bzw Delta^1/alpha #/ page 215
gamma_einsetzen <- Delta^(1/alpha) * sigma


Y_1 <- NULL
Y_2 <- NULL
N <- 100
i <- 0 



## loop 
while (length(Y_2) < N) {
  # generating random numbers
  U <- runif(1)
  V <- stabledist::rstable(1, alpha = alpha, beta = 1, gamma = gamma, delta = 0)
  
  # acceptance-rejection decision
  if(U <= exp(-b*(V + c))){
    Y_2 <- c(Y_2,  ( V - (delta_a * gamma(1 - alpha) * a * b^(alpha *-1))))
  } 
  # loop things
  i <- i + 1
  print(paste0('Y_2 = ', length(Y_2), '; ', 'tries = ', i))
}


# as a function 

tsd_algorithm_2 <- function(alpha, delta_a, a = 1, b = 1, c, N = 100){
  
  # c = constant; acts as a truncation
  # alpha parameter of the stable distribution 
  # a = skewness parameter delta_+ in LM 215: -> 1; Subordinator; absolute skewness -> beta in rstable()
  # b 
  # delta_a 
  # delta -> 0; location
  # gamma -> sigma^beta in LV and there delta_+ equals -> a in the simulation paper
  
  sigma_alpha <- (1/alpha) * gamma(1- alpha) * a * cos((alpha*pi)/2)
  gamma <- sigma_alpha^(1/alpha)
  
  # loop parameters
  Y_1 <- NULL
  Y_2 <- NULL
  i <- 0 
  
  while (length(Y_2) < N) {
    # generating random numbers
    U <- runif(1)
    V <- stabledist::rstable(1, alpha = alpha, beta = 1, gamma = gamma, delta = 0)
    
    # acceptance-rejection decision
    if(U <= exp(-b*(V + c))){
      Y_2 <- c(Y_2,  ( V - (delta_a * gamma(1 - alpha) * a * b^(alpha *-1))))
    } 
    # loop things
    i <- i + 1
    
    if (length(Y_2) %% 1000 == 0) {
      print(paste0('Y_2 = ', length(Y_2), '; ', 'tries = ', i))
    }
    
  }
  
  return(Y_2)
}

# possible to do it vector vise? 
# parallel
sim_v <- tsd_algorithm_2(alpha = 1.5, delta_a = 1, c = 0.06, N = 30000)
hist(sim_v)

#################################
#---- serial representation ----#
#################################
# parameters
k = 2 # k; number of series?
alpha <- 1.5
Delta <- 1
a <- 1
b <- 1
lambda_1 <- 1
lambda_2 <- 1

# random variables 
#pois_comp <- rpois(n = k, lambda = 1) # is standard 1? # arrival times of a poisson process not poisson itself
#arrival_pois <- cumsum(rexp(k)) # cumsum
arrival_gamma <- rgamma(k, shape = 1) # shape = 1, equals standard posssion process?


U <- runif(k)
E_1 <- rexp(k)
E_4 <- rexp(k) 
E_2 <- rexp(k, rate = b *lambda_1)
E_3 <- rgamma(k, shape = lambda_1, scale = (b * lambda_2)^(-1))

# 
gamma_delta <- (Delta * a / alpha)^(1/alpha) * VGAM::zeta(1/alpha) - Delta * gamma(1 - alpha) * a * b^(alpha - 1)

# computation of 5.2
sum( 
  # min of first part 
  pmin( 
((alpha* arrival_gamma)/ Delta* a)^(-1 / alpha), 
(E_1 * U^(1/alpha) /b)) - ((alpha* 1:k) / Delta * a)^(-1 / alpha))
  




#----- looop serial representation----
# parameters
N <- 100
Y <- rep(NA, N)
k <- 1e+06 # k; number of series
alpha <- 0.5
Delta <- 1
a <- 1
b <- 1
lambda_1 <- 1
lambda_2 <- 1


gamma_delta <- (Delta * a / alpha)^(1/alpha) * VGAM::zeta(1/alpha) - Delta * gamma(1 - alpha) * a * b^(alpha - 1)

for (i in seq_len(N)) {
  # random variables 
#  arrival_pois <- rexp(k) # 
  U <- runif(k)
  E_1 <- rexp(k)
  arrival_pois <- cumsum(rexp(k))
  
  # computation #of 5.2
  Y[i] <- sum( 
    # min of first part 
    pmin( 
      ((alpha* arrival_pois)/ Delta* a)^(-1 / alpha), 
      (E_1 * U^(1/alpha) /b)
      )
    )
# centering not needed in Subordinator 
# - ((alpha* 1:k) / Delta * a)^(-1 / alpha)
}

Y




#### function ####

serial_sim <- function(alpha, Delta, a = 1, b = 1, lambda_1 = 1, lambda_2 = 1, N = 100, k = 1e+04){
  
  # empty vector
  Y <- rep(NA, N)
  
  # simulation
  for (i in seq_len(N)) {
      # random variables 
      U <- runif(k)
      E_1 <- rexp(k)
      # arrival times of a Poisson process
      arrival_pois <- cumsum(rexp(k))
    
    # computation #of 5.2
      Y[i] <- sum(
        # min of first part 
        pmin( 
          ((alpha* arrival_pois)/ Delta* a)^(-1 / alpha), 
          (E_1 * U^(1/alpha) /b)
        )
    )
    
  }
  
  # Correction
  
  if(1 == abs(Delta)){
    gamma_delta <- (Delta * a / alpha)^(1/alpha) * VGAM::zeta(1/alpha) - Delta * gamma(1 - alpha) * a * b^(alpha - 1)
    
    centering <- ((alpha* 1:k) / Delta * a)^(-1 / alpha)
    
    # correction
    Y <- Y - centering + gamma_delta
  }
  
  return(Y)
}


try_pos <- serial_sim(alpha = 0.9, Delta = 0.5)

hist(try_pos)


# To do 
# -> arrival times of Poisson process ::: done 
# -> 2 versions! subordinator not subordinator version
# -> + gamma delta needed

# subordinator -> without centering part and without gamma (delta) (which is on the left side) 








#### not needed
# E_4 <- rexp(k) 
# E_2 <- rexp(k, rate = b *lambda_1)
# E_3 <- rgamma(k, shape = lambda_1, scale = (b * lambda_2)^(-1))
