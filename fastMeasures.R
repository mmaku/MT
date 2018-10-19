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
                     iterations = 250,  # Numer of graphs simulated to calculated FDR
                     epsilon = 10e-4,
                     testLocalFDR = TRUE)
{
   methods <- c("gLASSO", 
                 "banerjee.gLASSO",
                 "holm.gSLOPE",
                 "BH.gSLOPE")
    zeros <- c(0,0,0,0)
    results <- data.frame(procedure = methods, FDR = zeros, localFDR = zeros, Power = zeros, SP = zeros)
    
    # gLASSO parameters
    gLassoLambda <- lambdaGLASSO(p = p, n = n, alpha = alpha)
    banerjeeLambda <- lambdaBanerjee(p = p, n = n, alpha = alpha)
    
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
        
        omegaHat$GL     <- glasso(s = generatedData$sigmahat, rho = gLassoLambda, thr = epsilon,
                                  penalize.diagonal = penalizeDiagonal)$wi    
        omegaHat$GLB    <- glasso(s = generatedData$sigmahat, rho = banerjeeLambda, thr = epsilon,
                                  penalize.diagonal = penalizeDiagonal)$wi    
        omegaHat$GSH    <- gslopeADMM(sampleCovariance = generatedData$sigmahat, lambda = holmlambda,
                                      absoluteEpsilon = epsilon)$precisionMatrix
        omegaHat$GSBH   <- gslopeADMM(sampleCovariance = generatedData$sigmahat, lambda = BHlambda,
                                      absoluteEpsilon = epsilon)$precisionMatrix
        results$FDR         <- results$FDR + sapply(omegaHat, function(x) FDP(x, adjacent))
        if(testLocalFDR)
        {
            results$localFDR    <- results$FDR + sapply(omegaHat, function(x) localFDP(x, adjacent))    
        } else
        {
            results$localFDR    <- NA
        }
        results$Power       <- results$FDR + sapply(omegaHat, function(x) SN(x, adjacent))
        results$SP          <- results$FDR + sapply(omegaHat, function(x) SP(x, adjacent))
    }
    
    results[,-1] <- results[,-1]/iterations
    
    return(results)
}