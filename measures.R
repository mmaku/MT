# Written by Micha� Makowski

# install.packages("glasso")
# install.packages("huge")
# install.packages("R.utils")

require(glasso, quietly = TRUE)  # graphical lasso
require(huge, quietly = TRUE)    # graphical models
require(R.utils, quietly = TRUE) # doCall

source("auxilaryFunctions.R")
source("admmGSLOPE.R")
source("errorTypes.R")


measures <- function(n = 150, 
                     p = 200, 
                     graphType = "cluster",
                     graphParameters = NULL, 
                     alpha = 0.05, 
                     penalizeDiagonal = FALSE, 
                     additionalMethods = NULL,
                     iterations = 250,  # Numer of graphs simulated to calculated FDR
                     epsilon = 10e-4, 
                     verbose = TRUE)
{
    if(verbose) 
    {
        cat("Starting FDR, localFDR, Sensitivity & Specificity simulations\nn = ", 
            n, "\np = ", p, "\ngraph structure = ", graphType, 
            "\nsimulations number = ", iterations)
        
        progressBar <- txtProgressBar(min = 0, max = iterations, style = 3)
        setTxtProgressBar(progressBar, 0)
    }

    methods <- c("gLASSO", 
                 "banerjee.gLASSO",
                 "holm.gSLOPE",
                 "BH.gSLOPE", names(additionalMethods))
    zeros <- rep_len(0, length(methods))
    results <- data.frame(procedure = methods, FDR = zeros, localFDR = zeros, Power = zeros, SP = zeros)
    
    # gLASSO parameters
    gLassoLambda <- lambdaSelector(input = p, n = n, alpha = alpha, method = "glasso", verbose = FALSE)
    banerjeeLambda <- lambdaSelector(input = p, n = n, alpha = alpha, method = "banerjee", verbose = FALSE)
    
    # gSLOPE parameters
    holmlambda <- lambdaSelector(input = p, n = n, alpha = alpha, method = "holm", verbose = FALSE)
    BHlambda <- lambdaSelector(input = p, n = n, alpha = alpha, method = "BH", verbose = FALSE) 
    
    for(i in seq_len(iterations))
    {
        generatedData <- doCall("huge.generator", 
                                n = n, d = p, graph = graphType, verbose = FALSE, vis = FALSE,
                                args = graphParameters)
        
        adjacent <- properAdjacent(generatedData$theta)
        
        proc <- 1
        
        for(m in methods)
        {
            if(m == "gLASSO")
            {
                omegaHat <- glasso(s = generatedData$sigmahat, rho = gLassoLambda, thr = epsilon,
                                   penalize.diagonal = penalizeDiagonal)$wi    
            } else if(m == "banerjee.gLASSO")
            {
                omegaHat <- glasso(s = generatedData$sigmahat, rho = banerjeeLambda, thr = epsilon,
                                   penalize.diagonal = penalizeDiagonal)$wi    
            } else if(m == "holm.gSLOPE")
                
            {
                omegaHat <- gslopeADMM(sampleCovariance = generatedData$sigmahat, lambda = holmlambda,
                                       penalizeDiagonal = penalizeDiagonal,
                                       absoluteEpsilon = epsilon,
                                       verbose = FALSE)$precisionMatrix
            } else if(m == "BH.gSLOPE")
            {
                omegaHat <- gslopeADMM(sampleCovariance = generatedData$sigmahat, lambda = BHlambda,
                                       penalizeDiagonal = penalizeDiagonal, 
                                       absoluteEpsilon = epsilon, 
                                       verbose = FALSE)$precisionMatrix
            } else
            {
                parameters <- additionalMethods[[m]]
                additionalEstimation <- doCall("huge", x = generatedData$data, verbose = FALSE, args = parameters)
                
                omegaHat <- doCall("huge.select", est = additionalEstimation, verbose = FALSE, 
                                   args = parameters)$refit
                
                omegaHat <- properAdjacent(omegaHat)   
            }
        
            # omegaHat <- omegaHat*(omegaHat>0)
            
            results[proc,"FDR"] <- results[proc,"FDR"] + FDP(omegaHat, adjacent)
            results[proc,"localFDR"] <- results[proc,"localFDR"] + localFDP(omegaHat, adjacent)
            results[proc, "Power"] <- results[proc,"Power"] + SN(omegaHat, adjacent)
            results[proc, "SP"] <- results[proc,"SP"] + SP(omegaHat, adjacent)
            
            proc = proc + 1
        }

        if(verbose)
            setTxtProgressBar(progressBar, i)
    }
    
    results[,-1] <- results[,-1]/iterations
    
    if(verbose)
        close(progressBar)
    
    return(results)
}