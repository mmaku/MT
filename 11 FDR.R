# Written by Micha³ Makowski

# install.packages("glasso")
# install.packages("huge")

require(glasso)
require(huge)

source("01 auxilaryFunctions.R")
source("07 admmGSLOPE.R")
source("10 simpleMeasures.R")

measures <- function(n=150, 
                     p=200, 
                     graph="cluster",
                     alpha = .1, 
                     penalizeDiagonal = FALSE, 
                     epsilon = 10e-4, 
                     simulationsNumber = 1000,  # Numer of graphs simulated to calculated FDR
                     verbose = TRUE)
{
    if(verbose) 
    {
        cat("Starting FDR, Sensitivity & Specificity simulations\nn = ", 
            n, "\np = ", p, "\ngraph structure = ", graph, 
            "\nsimulations number = ", simulationsNumber, "\n")
        
        progressBar <- txtProgressBar(min = 0, max = simulationsNumber, style = 3)
        setTxtProgressBar(progressBar, 0)
    }
    
    gLASSO <- list(FDR = NULL,
                   SN = NULL,
                   SP = NULL)
    
    BSgSLOPE <- list(FDR = NULL,
                     SN = NULL,
                     SP = NULL)
    
    BHgSLOPE <- list(FDR = NULL,
                     SN = NULL,
                     SP = NULL)
    
    # gLASSO parameters
    banerjeeLambda <- lambdaSelector(p, n, alpha, method = "banerjee", verbose = FALSE)
        
    # gSLOPE parameters
    BSlambda <- lambdaSelector(p, n, alpha, method = "BS", verbose = FALSE) 
    BHlambda <- lambdaSelector(p, n, alpha, method = "BH", verbose = FALSE) 
    
    for(i in 1:simulationsNumber)
    {
        graphHUGE <- huge.generator(n, p, graph, verbose = FALSE) 
        
        omegaHATgLASSO <- glasso(graphHUGE$sigmahat, rho = banerjeeLambda, thr = epsilon,
                                 penalize.diagonal = penalizeDiagonal)$wi

        # omegaHATgLASSO <- glassoADMM(graphHUGE$sigmahat, lambda = banerjeeLambda, 
        #                              penalizeDiagonal = penalizeDiagonal, 
        #                              truncate = TRUE, absoluteEpsilon = epsilon, verbose = FALSE)$precisionMatrix
        
        BSomegaHATgSLOPE <- gslopeADMM(graphHUGE$sigmahat, lambda = BSlambda,
                                       penalizeDiagonal = penalizeDiagonal, 
                                       truncate = TRUE, absoluteEpsilon = epsilon, verbose = FALSE)$precisionMatrix

        BHomegaHATgSLOPE <- gslopeADMM(graphHUGE$sigmahat, lambda = BHlambda,
                                       penalizeDiagonal = penalizeDiagonal, 
                                       truncate = TRUE, absoluteEpsilon = epsilon, verbose = FALSE)$precisionMatrix
        
        gLASSO$FDR <- c(gLASSO$FDR, FDP(omegaHATgLASSO, graphHUGE$omega))
        gLASSO$SN  <- c(gLASSO$SN, SN(omegaHATgLASSO, graphHUGE$omega))
        gLASSO$SP  <- c(gLASSO$SP, SP(omegaHATgLASSO, graphHUGE$omega))
        
        BSgSLOPE$FDR <- c(BSgSLOPE$FDR, FDP(BSomegaHATgSLOPE, graphHUGE$omega))
        BSgSLOPE$SN  <- c(BSgSLOPE$SN, SN(BSomegaHATgSLOPE, graphHUGE$omega))
        BSgSLOPE$SP  <- c(BSgSLOPE$SP, SP(BSomegaHATgSLOPE, graphHUGE$omega))
        
        BHgSLOPE$FDR <- c(BHgSLOPE$FDR, FDP(BHomegaHATgSLOPE, graphHUGE$omega))
        BHgSLOPE$SN  <- c(BHgSLOPE$SN, SN(BHomegaHATgSLOPE, graphHUGE$omega))
        BHgSLOPE$SP  <- c(BHgSLOPE$SP, SP(BHomegaHATgSLOPE, graphHUGE$omega))
        
        if(verbose)
            setTxtProgressBar(progressBar, i)
    }
    
    if(verbose)
        close(progressBar)
    
    return(list(gLASSO = lapply(gLASSO, function(x) sum(x)/simulationsNumber),
                BSgSLOPE = lapply(BSgSLOPE, function(x) sum(x)/simulationsNumber),
                BHgSLOPE = lapply(BHgSLOPE, function(x) sum(x)/simulationsNumber)))
}
