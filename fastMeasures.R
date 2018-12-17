# Written by Micha³ Makowski

# install.packages("glasso")
# install.packages("huge")
# install.packages("R.utils")

require(glasso, quietly = TRUE)  # graphical lasso
require(huge, quietly = TRUE)    # graphical models
require(R.utils, quietly = TRUE) # doCall

source("auxilaryFunctions.R")
source("fastadmmGSLOPE.R")
source("errorTypes.R")

measures <- function(n = 150, 
                     p = 200, 
                     graphType = "cluster",
                     graphParameters = NULL, 
                     alpha = 0.05, 
                     penalizeDiagonal = FALSE, 
                     scaled = TRUE,
                     iterations = 250,  # Numer of graphs simulated to calculated FDR
                     epsilon = 10e-4,
                     testLocalFDR = TRUE)
{
    methods <- c("gLASSO",
                 "banerjee.gLASSO",
                 "holm.gSLOPE",
                 "BH.gSLOPE")
    zeros <- rep_len(0, length.out = length(methods))
    results <- data.frame(procedure = methods, FDR = zeros, localFDR = zeros, 
                          Power = zeros, SP = zeros, algIter = zeros) # algIter for benchmarking
    
    # Fixing bad implementation of huge
    if(!is.null(graphParameters$u))
        graphParameters$u <- graphParameters$u - 0.1
    
    # gLASSO parameters
    banerjeeLambda <- lambdaBanerjee(p = p, n = n, alpha = alpha)
    gLassoLambda <- lambdaGLASSO(p = p, n = n, alpha = alpha)
    
    # gSLOPE parameters
    holmlambda <- lambdaHolm(p = p, n = n, alpha = alpha)
    BHlambda <- lambdaBH(p = p, n = n, alpha = alpha)
    
    zeros <- matrix(0, nrow = p ,ncol = p)
    omegaHat <- list(GL = zeros, GLB = zeros, GSH = zeros, GSBH = zeros)
    
    
    
    for(i in seq_len(iterations))
    {
        generatedData <- doCall("huge.generator", 
                                n = n, d = p, graph = graphType, verbose = FALSE, vis = FALSE,
                                args = graphParameters)
        
        adjacent        <- properAdjacent(generatedData$theta)
        
        if(scaled)
        {
            initialMatrix <- generatedData$sigmahat # correlation matrix
        } else
        {
            initialMatrix <- cov(generatedData$data)
        }
        
        GL    <- glasso(s = initialMatrix, rho = gLassoLambda, thr = epsilon,
                         penalize.diagonal = penalizeDiagonal)
        omegaHat$GL <- GL$wi
        
        GLB    <- glasso(s = initialMatrix, rho = banerjeeLambda, thr = epsilon,
                         penalize.diagonal = penalizeDiagonal)
        omegaHat$GLB <- GLB$wi
        
        GSH    <- gslopeADMM(sampleCovariance = initialMatrix, lambda = holmlambda,
                             Y = omegaHat$GLB,
                             absoluteEpsilon = epsilon)
        omegaHat$GSH <- GSH$precisionMatrix
        
        GSBH   <- gslopeADMM(sampleCovariance = initialMatrix, lambda = BHlambda,
                             Y = omegaHat$GSH,
                             absoluteEpsilon = epsilon)
        omegaHat$GSBH <- GSBH$precisionMatrix
        
        results$FDR         <- results$FDR + sapply(omegaHat, function(x) FDP(x, adjacent))
        if(testLocalFDR) # localFDR uses functions from igraph that are expansive computationally
        {
            if((generatedData$graph.type == "cluster") & isTRUE(graphParameters$prob == 1))
            {
                results$localFDR <- results$FDR
            } else
            {
                results$localFDR <- results$localFDR + sapply(omegaHat, function(x) localFDP(x, adjacent))       
            }
        } else
        {
            results$localFDR    <- NA
        }
        results$Power       <- results$Power + sapply(omegaHat, function(x) SN(x, adjacent))
        results$SP          <- results$SP + sapply(omegaHat, function(x) SP(x, adjacent))
        
        results$algIter  <- results$algIter + c(GL$niter, GLB$niter, GSH$iterations, GSBH$iterations)
    }
    
    results[,-1] <- results[,-1]/iterations
    
    return(results)
}