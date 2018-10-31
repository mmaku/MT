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

createRocMatrix <- function(nVec = 150, 
                                   pVec = 200, 
                                   graphTypeVec = "cluster",
                                   alphaVec = 0.05,
                                   iterationsVec = 100)
{
    output <- expand.grid(nVec, pVec, graphTypeVec, alphaVec, penalizeDiagonalVec, scaledVec, iterationsVec,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    colnames(output) <- c("n", "p", "graphType", "alpha", "penalizeDiagonal", "scaled", "iterations")
    
    return(output)
}


graphROC <- function(n = 150, 
                     p = 200, 
                     graphType = "cluster",
                     graphParameters = NULL, 
                     alpha = 0.05, 
                     alphaMultiSeq = seq(from=1/2, to=3/2, length.out=10),
                     penalizeDiagonal = FALSE, 
                     scaled = TRUE,
                     iterations = 250,
                     epsilon = 10e-4,
                     verbose = TRUE)
{
    if(verbose) 
    {
        cat("Starting ROC simulations\nn = ", 
            n, "\np = ", p, "\ngraph structure = ", graphType, 
            "\nsimulations number = ", iterations)
        
        progressBar <- txtProgressBar(min = 0, max = iterations, style = 3)
        setTxtProgressBar(progressBar, 0)
    }
    # Fixing bad implementation of huge
    if(!is.null(graphParameters$u))
        graphParameters$u <- graphParameters$u - 0.1
    
    # gLASSO parameters
    banerjeeLambda <- lambdaBanerjee(p = p, n = n, alpha = alpha)
    
    # gSLOPE parameters
    holmlambda <- lambdaHolm(p = p, n = n, alpha = alpha)
    BHlambda <- lambdaBH(p = p, n = n, alpha = alpha)
    
    # omegaHat list
    zeros <- matrix(0, nrow = p ,ncol = p)
    omegaHat <- list(GLB = zeros, GSH = zeros, GSBH = zeros)
    
    # results array
    methods <- c("banerjee.gLASSO",
                 "holm.gSLOPE",
                 "BH.gSLOPE")
    
    measures <- c("multiplier", "SN", "SP")
    
    results <- array(0, 
                     dim = c(3, length(alphaMultiSeq), 3), 
                     dimnames = list(measures, NULL, methods))
    
    results["multiplier",,] <- alphaMultiSeq
    
    for(i in seq_len(iterations))
    {
        generatedData <- doCall("huge.generator", 
                                n = n, d = p, graph = graphType, verbose = FALSE,
                                args = graphParameters)
        
        adjacent <- properAdjacent(generatedData$theta)
        
        if(scaled)
        {
            initialMatrix <- generatedData$sigmahat # correlation matrix
        } else
        {
            initialMatrix <- cov(generatedData$data)
        }
        
        for(multiIter in seq_along(alphaMultiSeq))
        {
            GLB    <- glasso(s = initialMatrix, rho = banerjeeLambda*alphaMultiSeq[multiIter], 
                             thr = epsilon,
                             penalize.diagonal = penalizeDiagonal)
            omegaHat$GLB <- GLB$wi
            
            GSH    <- gslopeADMM(sampleCovariance = initialMatrix, 
                                 lambda = holmlambda*alphaMultiSeq[multiIter],
                                 Y = omegaHat$GLB,
                                 absoluteEpsilon = epsilon)
            omegaHat$GSH <- GSH$precisionMatrix
            
            GSBH   <- gslopeADMM(sampleCovariance = initialMatrix, 
                                 lambda = BHlambda*alphaMultiSeq[multiIter],
                                 Y = omegaHat$GSH,
                                 absoluteEpsilon = epsilon)
            omegaHat$GSBH <- GSBH$precisionMatrix
            
            results["SN", multiIter,] <- results["SN", multiIter,] + 
                sapply(omegaHat, function(x) SN(x, adjacent))
            results["SP", multiIter,] <- results["SP", multiIter,] + 
                sapply(omegaHat, function(x) SP(x, adjacent))
        }
        
        if(verbose)
            setTxtProgressBar(progressBar, i)
    }
    
    results[-1,,] <- results[-1,,]/iterations
    
    if(verbose)
        close(progressBar)
    
    return(results)
}
