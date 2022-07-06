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
# a > 0 muss größer null sein

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
#2  functions -> subordinator
# a > 0 muss größer null sein, Scaling property?




sub_serial_sim <- function(alpha, Delta = 1, a, b, N = 100, k = 1e02){
  
  if(alpha <= 0 | alpha >= 1) 
    stop("'alpha' must be in ]0,1[")
  # empty vector
  Y <- rep(NA, N)
    # simulation
    for (i in seq_len(N)){
      # random variables 
      U <- runif(k)
      E_1 <- rexp(k)
      # arrival times of a Poisson process
      arrival_pois <- cumsum(rexp(k))
      
      # computation #of 5.2
      Y[i] <- sum(
        # min of first part 
        pmin( 
          ((alpha* arrival_pois)/ Delta*a)^(-1 / alpha), 
          ((E_1 * U^(1/alpha)) /b)
        )
      )
    }
  return(Y)
} 

try_1 <- sub_serial_sim(alpha = 0.9, Delta = 1, a = 1, b = 1, N = 1e6)

hist(try_1, breaks = 100)



#### function ####
serial_sim_fun <- function(alpha, Delta, a = 1, b = 1, N = 100, k = 1e+04){
  if(alpha >= 2 | alpha <= 1) 
    stop("'alpha' must be in ]1,2[")
  
  # empty vector
  Y <- rep(NA, N)

  # correction terms 
  gamma_delta <- ((Delta * a / alpha)^(1/alpha)) * VGAM::zeta((1/alpha)) - (Delta * gamma(1 - alpha) * a * b^(alpha - 1))
  centering <- ((alpha* 1:k) / Delta * a)^(-1 / alpha)
    
  # simulation
  for (i in seq_len(N)){
    # random variables 
    U <- runif(k)
    E_1 <- rexp(k)
    # arrival times of a Poisson process
    arrival_pois <- cumsum(rexp(k))
      
      # computation #of 5.2
    Y[i] <- sum(
      # min of first part 
      (pmin( 
        ((alpha* arrival_pois)/ Delta*a)^(-1 / alpha), 
        ((E_1 * U^(1/alpha)) /b)
      )
        ###### Correction ######
        - centering) 
        
      )
  }
  
  Y <- Y + gamma_delta
  
  return(Y)
}

try_serial <- serial_sim_fun(alpha = 1.0, Delta = 1, a = 1, b = 1, N = 1e1)

hist(try_serial, breaks = 100)



serial_sim <- function(alpha, Delta = 1, a_p = 1, b_p = 1, a_m = 1, b_m = 1, N = 1e02, k = 1e04){
  
  if (alpha > 2 | alpha <= 1) {
    
  }
  
  x_p <- serial_sim_fun(alpha = alpha, Delta = Delta, a = a_p, b = b_p, N = N)
  
  x_m <- serial_sim_fun(alpha = alpha, Delta = Delta, a = a_m, b = b_m, N = N)
  
  x = x_p - x_m
  
  return(x)
  
}


try <- serial_sim(alpha = 1.2, N = 1e04)

hist(try)




# which Parameter does what?
# Delta -> stepsize -> 1 standard
# alpha, alpha or beta parameter of alpha/beta stable distribution