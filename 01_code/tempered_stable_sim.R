## ---- packages ----
source(here::here('00_packages/packages.R'))

## ---- simulation ----
# On simulation of tempered stable random variates


## --- Algorithm 0 ----


b <- 0.1
alpha <- 0.2
# beta -> 1; Subordinator; absolut skewness; 
# delta -> 0; location
# gamma -> 


U <- runif(1)

# Levy matters IV -> page 215
V <- stabledist::rstable(1, alpha = alpha, beta = 1, gamma = 1, delta = 0)


if(U <= exp(-b*V)) Y <- V


