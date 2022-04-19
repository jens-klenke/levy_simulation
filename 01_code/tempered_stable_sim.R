## ---- packages ----
source(here::here('00_packages/packages.R'))

## ---- simulation ----
# On simulation of tempered stable random variates


## --- Algorithm 0 ----

U <- runif(1)

alpha <- 0.2
# beta -> 1; Subordinator; absolut skewness; 

# delta -> 0; location
# Levy matters IV -> page 215
S <- stabledist::rstable(1, alpha = alpha, gamma = , beta = 1, delta = 0)



