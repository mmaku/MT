# Written by Micha³ Makowski

# install.packages("SLOPE")

source("06 OWL1prox.R")

require(SLOPE)

# C - sample covariance
# mu - Augmented Lagrangian parameter
# lambda - Regularization parameter controlling sparsity
# maxIter - maximum number of iterations
# absoluteEpsilon - used in residual stopping criterium
# truncate - should entries below absoluteEpsilon be equal to zero?
gslopeADMM <- function(sampleCovariance, 
                       mu = 1, 
                       lambda = NULL, 
                       penalizeDiagonal = TRUE,
                       maxIter = 1e5, 
                       absoluteEpsilon = 1e-4, 
                       relativeEpsilon = 1e-4, 
                       truncate = F)
{
    lambda <- sort(lambda, decreasing = T)
    
    entriesNumber <- sum(1:(ncol(sampleCovariance)-!penalizeDiagonal))
    
    if(length(lambda) < entriesNumber) 
        lambda <- c(lambda, rep(0, times = entriesNumber - length(lambda)))
    
    if(length(lambda) > entriesNumber) 
        lambda <- lambda[1:entriesNumber]   
    
    
    Z <- sampleCovariance*0 # Initialize Lagragian to be nothing (seems to work well)
    Y <- Z 
    X <- diag(nrow = nrow(sampleCovariance))

    
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

        if(primalResidual < primalEpsilon & dualResidual < dualEpsilon) 
            break
        
        if(n %% 10000 == 0)
            print(paste(n, "iterations done..."), quote = F)
    }
    
    if(truncate)
        X[X < absoluteEpsilon] <- 0
    
    print("ADMM gSLOPE done.")
    
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
