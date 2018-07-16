# Written by Micha³ Makowski

# install.packages("glasso")
# install.packages("huge")
# install.packages("R.utils")

require(glasso)  # graphical lasso
require(huge)    # graphical models
require(R.utils) # doCall

source("01 auxilaryFunctions.R")
source("07 admmGSLOPE.R")
source("10 simpleMeasures.R")

measures <- function(n = 150, 
                     p = 200, 
                     graphType = "cluster",
                     alpha = 0.1, 
                     penalizeDiagonal = FALSE, 
                     partialMeasure = TRUE,
                     additionalMethods = NULL,
                     iterations = 1000,  # Numer of graphs simulated to calculated FDR
                     epsilon = 10e-4, 
                     verbose = TRUE)
{
    if(verbose) 
    {
        cat("Starting FDR, Sensitivity & Specificity simulations\nn = ", 
            n, "\np = ", p, "\ngraph structure = ", graphType, 
            "\nsimulations number = ", iterations, 
            "\nmeasures calculated on whole matrix = ", c("No.\n", "Yes.\n")[1+partialMeasure])
        
        progressBar <- txtProgressBar(min = 0, max = iterations, style = 3)
        setTxtProgressBar(progressBar, 0)
    }

    methods <- c("gLASSO", "BSgSLOPE", "BHgSLOPE", names(additionalMethods))
    zeros <- rep_len(0, length(methods))
    results <- data.frame(FDR = zeros, SN = zeros, SP = zeros, row.names = methods)
   
    # gLASSO parameters
    banerjeeLambda <- lambdaSelector(input = p, n = n, alpha = alpha, method = "banerjee", verbose = FALSE)
        
    # gSLOPE parameters
    BSlambda <- lambdaSelector(input = p, n = n, alpha = alpha, method = "BS", verbose = FALSE) 
    BHlambda <- lambdaSelector(input = p, n = n, alpha = alpha, method = "BH", verbose = FALSE) 
    
    for(i in 1:iterations)
    {
        generatedData <- huge.generator(n, d = p, graph = graphType, verbose = FALSE) 
        adjacent <- properAdjacent(generatedData$theta)
        
        for(m in methods)
        {
            if(m == "gLASSO")
            {
                omegaHat <- glasso(s = generatedData$sigmahat, rho = banerjeeLambda, thr = epsilon,
                                   penalize.diagonal = penalizeDiagonal)$wi    
            } else if(m == "BSgSLOPE")
            {
                omegaHat <- gslopeADMM(sampleCovariance = generatedData$sigmahat, lambda = BSlambda,
                                       penalizeDiagonal = penalizeDiagonal, 
                                       truncate = TRUE, absoluteEpsilon = epsilon, 
                                       verbose = FALSE)$precisionMatrix
            } else if(m == "BHgSLOPE")
            {
                omegaHat <- gslopeADMM(sampleCovariance = generatedData$sigmahat, lambda = BHlambda,
                                       penalizeDiagonal = penalizeDiagonal, 
                                       truncate = TRUE, absoluteEpsilon = epsilon, 
                                       verbose = FALSE)$precisionMatrix
            } else
            {
                parameters <- additionalMethods[[m]]
                additionalEstimation <- doCall("huge", x = generatedData$data, verbose = FALSE, args = parameters)
                
                omegaHat <- doCall("huge.select", est = additionalEstimation, verbose = FALSE, 
                                   args = parameters)$refit
                
                omegaHat <- properAdjacent(omegaHat)   
            }
            
            results[m,"FDR"] <- results[m,"FDR"] + FDP(omegaHat, adjacent, partialMeasure)
            results[m,"SN"] <- results[m,"SN"] + SN(omegaHat, adjacent, partialMeasure)
            results[m,"SP"] <- results[m,"SP"] + SP(omegaHat, adjacent, partialMeasure)
        }

        if(verbose)
            setTxtProgressBar(progressBar, i)
    }
    
    if(verbose)
        close(progressBar)
    
    return(results/iterations)
}