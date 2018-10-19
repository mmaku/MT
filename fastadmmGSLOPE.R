# Written by Micha³ Makowski

# TODO
# Checks
# Lamda default value

# install.packages("SLOPE")
require(SLOPE)

source("OWL1prox.R")

# sampleCovariance - sample covariance
# mu - Augmented Lagrangian parameter
# lambda - sequence of lambda regularizators (ordered L1 norm)
# maxIter - maximum number of iterations
# absoluteEpsilon - used in residual stopping criterium
# verbose - console output

gslopeADMM <- function(sampleCovariance, 
                       lambda = NULL,
                       maxIter = 1e5, 
                       absoluteEpsilon = 1e-4)
{
    p <- ncol(sampleCovariance)
    
    # Sequence length
    
    entriesNumber <- sum(1:(p-1))
    
    # Lambda preparation
    
    lambda <- lambda[seq(from = p+1, to = length(lambda), by = 2)]

    # Initialization
    
    Z <- sampleCovariance*0 # Initialize Lagragian to be nothing (seems to work well)
    Y <- Z 
    X <- diag(nrow = p)

    # ADMM algotithm
    
    for(n in 1:maxIter)
    {
        # Solve sub-problem to solve X
        Ctilde <- Y-Z-sampleCovariance
        Ceigen <- eigen(Ctilde, symmetric = TRUE)
        CeigenVal <- Ceigen$val
        CeigenVec <- Ceigen$vec
        Fmu <- 1/2*diag(CeigenVal+sqrt(CeigenVal*CeigenVal+4))
        X <- CeigenVec%*%Fmu%*%t(CeigenVec)
        
        # Solve sub-problem to solve Y
        Yold <- Y 
        # Y <- softThresholding(X+Z, lambda/mu) 
        Y <- matrixOWL1prox(X+Z, lambda, FALSE) 
        
        # Update the Lagrangian
        Z <- Z + (X-Y)
        
        # Residuals
        primalResidual <- norm(X-Y, type = "F")
        dualResidual   <- norm((Y-Yold), type = "F")
        
        # Stopping criteria
        primalEpsilon <- absoluteEpsilon # + relativeEpsilon*max(l2norm(X), l2norm(Y))
        # dualEpsilon   <- absoluteEpsilon # + relativeEpsilon*l2norm(Z)

        if(primalResidual < primalEpsilon & dualResidual < primalEpsilon) 
            break
    }
    
    X[abs(X) < absoluteEpsilon] <- 0
    
    return(list(sampleCovariance = sampleCovariance,
                lambda = lambda, 
                precisionMatrix = X, 
                covarianceMatrix = solve(X), 
                residuals = c(primalResidual, dualResidual), 
                iterations = n, 
                epsilon = absoluteEpsilon))
}
