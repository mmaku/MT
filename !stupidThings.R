# Written by Micha³ Makowski

# install.packages("tictoc")
# install.packages("rbenchmark")
require(tictoc)
require(rbenchmark)

source("softThresholding.R")


tic("total")
tic("init")

gglasso <- function(x, lam, diag = F)
{
    glasso::glasso(x, rho = lam, penalize.diagonal = diag)
}

hglasso <- function(x, lam, scr = F)
{
    huge.glasso(graphHUGE$sigmahat, lambda = lam, scr = scr, verbose = FALSE)
}

huge.select(huge(graphHUGE$sigmahat, method = "ct")


benchmark(gglasso(graphHUGE$sigmahat, lam =  banerjeeLambda), 
          hglasso(graphHUGE$sigmahat, lam = banerjeeLambda), 
          replications = 10)


   