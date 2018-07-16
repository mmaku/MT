# Written by Micha³ Makowski

# TODO
# Checks
# Lamda default value

source("02 softThresholding.R")

# sampleCovariance - sample covariance
# mu - Augmented Lagrangian parameter
# lambda - Regularization parameter controlling sparsity
# alpha - if lambda is NULL, then alpha is used for sequence construction
# maxIter - maximum number of iterations
# absoluteEpsilon - used in residual stopping criterium
# truncate - should entries below absoluteEpsilon be equal to zero?
# verbose - console output

glassoADMM <- function(sampleCovariance, 
                       mu = 1, 
                       lambda = 1,
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
        cat("Starting ADMM gLASSO procedure...")
        progressBar <- txtProgressBar(min = 0, max = 1/absoluteEpsilon, style = 3)
        setTxtProgressBar(progressBar, 0)
    }
    
    # Checks
    
    if(!is.matrix(lambda) & length(lambda)!=1 & length(lambda)!=nrow(sampleCovariance))
    {
        stop("Wrong number of elements in lambda")
    }
    
    if(length(lambda)==1)
    {
        if(lambda==0)
        { 
            warning("With lambda=0, there may be convergence problems if the input matrix is not of full rank")
        }
    }
    
    # Lambda preparation
    
    lambda <- round(lambda, round(log10((1/absoluteEpsilon))))
    
    # Initialization
    
    Z <- sampleCovariance*0 # Initialize Lagragian to be nothing (seems to work well)
    Y <- Z 
    X <- diag(nrow = nrow(sampleCovariance), ncol = ncol(sampleCovariance))

    # ADMM algorithm
    
    for(n in 1:maxIter)
    {
        # Solve sub-problem to solve X
        Ctilde <- Y-Z-sampleCovariance/mu
        Ceigen <- eigen(Ctilde)
        CeigenVal <- Ceigen$val
        CeigenVec <- Ceigen$vec
        Fmu <- 1/2*diag(CeigenVal+sqrt(CeigenVal*CeigenVal+4/mu))
        X <- CeigenVec%*%Fmu%*%t(CeigenVec)
        
        # Solve sub-problem to solve Y
        Yold <- Y 
        Y <- softThresholding(X+Z, lambda/mu, penalizeDiagonal) 
        
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
