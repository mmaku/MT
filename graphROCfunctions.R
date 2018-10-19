# Written by Micha³ Makowski

# install.packages("glasso")
# install.packages("huge")
# install.packages("R.utils")

require(glasso, quietly = TRUE)  # graphical lasso
require(huge, quietly = TRUE)    # graphical models
require(R.utils, quietly = TRUE) # doCall

source("auxilaryFunctions.R")
source("admmGSLOPE.R")
source("errorTypes.R")


graphROC <- function(n = 150, 
                     p = 200, 
                     graphType = "cluster",
                     graphParameters = NULL, 
                     method = "banerjee",
                     alpha = 0.05, 
                     alphaMultiSeq = seq(from=1/2, to=3/2, length.out=10),
                     penalizeDiagonal = FALSE, 
                     iterations = 250,  # Numer of graphs simulated to calculated FDR
                     epsilon = 10e-4, 
                     verbose = TRUE)
{
    if(verbose) 
    {
        cat("Starting ROC simulations\nn = ", 
            n, "\np = ", p, "\ngraph structure = ", graphType, "\nalpha = ", alpha, 
            "\nsimulations number = ", iterations)
        
        progressBar <- txtProgressBar(min = 0, max = iterations, style = 3)
        setTxtProgressBar(progressBar, 0)
    }

    # parameters
    lambda <- lambdaSelector(input = p, n = n, alpha = alpha, method = method, verbose = FALSE)
    
    methods <- paste0("ROC", method)
    zeros <- rep_len(0, length(alphaMultiSeq))
    results <- data.frame(procedure = methods, multi = alphaMultiSeq, FDR = zeros, localFDR = zeros, SN = zeros, SP = zeros)

    
    for(i in seq_len(iterations))
    {
        generatedData <- doCall("huge.generator", 
                                n = n, d = p, graph = graphType, verbose = FALSE,
                                args = graphParameters)
        
        adjacent <- properAdjacent(generatedData$theta)
        
        proc <- 1
        
        for(multi in alphaMultiSeq)
        {
            if(method == "banerjee")
            {
                omegaHat <- glasso(s = generatedData$sigmahat, rho = lambda*multi, thr = epsilon,
                                   penalize.diagonal = penalizeDiagonal)$wi    
            } else
            {
                omegaHat <- gslopeADMM(sampleCovariance = generatedData$sigmahat, lambda = lambda*multi,
                                       penalizeDiagonal = penalizeDiagonal,
                                       absoluteEpsilon = epsilon,
                                       verbose = FALSE)$precisionMatrix
            }
            
            # omegaHat <- omegaHat*(omegaHat>0)
            
            results[proc,"FDR"] <- results[proc,"FDR"] + FDP(omegaHat, adjacent)
            results[proc,"localFDR"] <- results[proc,"localFDR"] + localFDP(omegaHat, adjacent)
            results[proc, "SN"] <- results[proc,"SN"] + SN(omegaHat, adjacent)
            results[proc, "SP"] <- results[proc,"SP"] + SP(omegaHat, adjacent)
            
            proc = proc + 1
        }
        
        if(verbose)
            setTxtProgressBar(progressBar, i)
    }
    
    results[,c(-1,-2)] <- results[,c(-1,-2)]/iterations
    
    if(verbose)
        close(progressBar)
    
    return(results)
}
