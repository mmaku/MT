# Written by Micha³ Makowski

# TODO
# Checks
# Lamda default value

# install.packages("SLOPE")
require(SLOPE)

source("06 OWL1prox.R")

# sampleCovariance - sample covariance
# mu - Augmented Lagrangian parameter
# lambda - sequence of lambda regularizators (ordered L1 norm)
# alpha - if lambda is NULL, then alpha is used for sequence construction
# maxIter - maximum number of iterations
# absoluteEpsilon - used in residual stopping criterium
# truncate - should entries below absoluteEpsilon be equal to zero?
# verbose - console output

gslopeADMM <- function(sampleCovariance, 
                       mu = 1, 
                       lambda = NULL, 
                       alpha = NULL,
                       penalizeDiagonal = TRUE,
                       maxIter = 1e5, 
                       absoluteEpsilon = 1e-4, 
                       relativeEpsilon = 1e-4, 
                       truncate = TRUE,
                       verbose = TRUE)
{
    # Console output
    
    if(verbose) 
    {
        cat("Starting ADMM gsLOPE procedure...")
        progressBar <- txtProgressBar(min = 0, max = 1/absoluteEpsilon, style = 3)
        setTxtProgressBar(progressBar, 0)
    }    
    
    p <- ncol(sampleCovariance)
    
    # Sequence length
    
    entriesNumber <- sum(1:(p-!penalizeDiagonal))
    
    # Lambda preparation
    
    lambda <- sort(lambda, decreasing = T)

    if(length(lambda) == p^2)
    {
        if(penalizeDiagonal)
        {
            lambda <- c(lambda[1:p], lambda[seq(from = p+1, to = length(lambda), by = 2)])
        } else
        {
            lambda <- lambda[seq(from = p+1, to = length(lambda), by = 2)]
        }
    } else if(length(lambda) < entriesNumber) 
    {
        lambda <- c(lambda, rep(0, times = entriesNumber - length(lambda)))
    } else if(length(lambda) > entriesNumber) 
    {
        lambda <- lambda[1:entriesNumber]   
    }
        
    # Initialization
    
    Z <- sampleCovariance*0 # Initialize Lagragian to be nothing (seems to work well)
    Y <- Z 
    X <- diag(nrow = p)

    # ADMM algotithm
    
    for(n in 1:maxIter)
    {
        # Solve sub-problem to solve X
        Ctilde <- Y-Z-sampleCovariance/mu
        Ceigen <- eigen(Ctilde, symmetric = TRUE)
        CeigenVal <- Ceigen$val
        CeigenVec <- Ceigen$vec
        Fmu <- 1/2*diag(CeigenVal+sqrt(CeigenVal*CeigenVal+4/mu))
        X <- CeigenVec%*%Fmu%*%t(CeigenVec)
        
        # Solve sub-problem to solve Y
        Yold <- Y 
        # Y <- softThresholding(X+Z, lambda/mu) 
        Y <- matrixOWL1prox(X+Z, lambda/mu, penalizeDiagonal) 
        
        # Update the Lagrangian
        Z <- Z + mu*(X-Y)
        
        # Residuals
        primalResidual <- norm(X-Y, type = "F")
        dualResidual   <- norm(mu*(Y-Yold), type = "F")
        
        # Stopping criteria
        primalEpsilon <- absoluteEpsilon # + relativeEpsilon*max(l2norm(X), l2norm(Y))
        dualEpsilon   <- absoluteEpsilon # + relativeEpsilon*l2norm(Z)

        if(verbose)
            setTxtProgressBar(progressBar, min(1/primalResidual, 1/dualResidual, 1/absoluteEpsilon))
        
        if(primalResidual < primalEpsilon & dualResidual < dualEpsilon) 
            break
    }
    
    if(truncate)
        X[X < absoluteEpsilon] <- 0
    
    if(verbose) 
        close(progressBar)
    
    return(list(sampleCovariance = sampleCovariance,
                lambda = lambda, 
                lagrangianParameter = mu,
                diagonalPenalization = penalizeDiagonal,
                precisionMatrix = X, 
                covarianceMatrix = solve(X), 
                residuals = c(primalResidual, dualResidual), 
                iterations = n, 
                truncated = truncate,
                epsilon = absoluteEpsilon))
}
