# Written by Micha³ Makowski

# install.packages("MASS")
# install.packages("glasso")
# install.packages("huge")

require(glasso)
require(huge)

source("01 auxilaryFunctions.R")
source("04 admmGLASSO.R")
source("07 admmGSLOPE.R")
source("10 contingencyMeasures.R")

FDR <- function(n=200, 
                d=200, 
                graph="cluster",
                alpha = .05, 
                penalizeDiagonal = FALSE, 
                epsilon = 10e-4, 
                simulationsNumber = 1000,  # Numer of graphs simulated to calculated FDR
                verbose = TRUE)
{
    if(verbose) 
    {
        cat("Starting FDR simulations\nn = ", n, "\nd = ", d, "\ngraph structure = ", graph, 
                    "\nsimulations number = ", simulationsNumber, "\n")
        progressBar <- txtProgressBar(min = 1, max = simulationsNumber,style = 3)
        setTxtProgressBar(progressBar, 1)
    }
    
    FDPgLASSO <- NULL
    FDPgSLOPE <- NULL
    
    SNgLASSO <- NULL
    SNgSLOPE <- NULL
    
    # gLASSO parameters
    m     <- d*(d-1)/2
    mBanerjee <- d^2   
    banerjeeLassoLambda <- qt(1-alpha/2/mBanerjee, df = n-2)/sqrt(n-2+qt(1-alpha/2/mBanerjee, df = n-2)^2)
    
    # gSLOPE parameters
    k = 1:m
    BHlambda <- qt(1-alpha*k/2/m, df = n-2)/sqrt(n-2+qt(1-alpha*k/2/m, df = n-2)^2)
    BHSlopeLambda <- c(rep(BHlambda[1], d), rep(BHlambda, each=2))
    
    for(i in 1:simulationsNumber)
    {
        graphHUGE <- huge.generator(n, d, graph, verbose = FALSE) 
        
        omegaHATgLASSO <- glasso(graphHUGE$sigmahat, rho = banerjeeLassoLambda, thr = epsilon,
                                 penalize.diagonal = penalizeDiagonal)$wi

        # omegaHATgLASSO <- glassoADMM(graphHUGE$sigmahat, lambda = banerjeeLassoLambda, 
        #                              penalizeDiagonal = penalizeDiagonal, 
        #                              truncate = TRUE, absoluteEpsilon = epsilon, verbose = FALSE)$precisionMatrix
        
        omegaHATgSLOPE <- gslopeADMM(graphHUGE$sigmahat, lambda = BHSlopeLambda, 
                                     penalizeDiagonal = penalizeDiagonal, 
                                     truncate = TRUE, absoluteEpsilon = epsilon, verbose = FALSE)$precisionMatrix
        
        FDPgLASSO <- c(FDPgLASSO, FDP(omegaHATgLASSO, graphHUGE$omega))
        FDPgSLOPE <- c(FDPgSLOPE, FDP(omegaHATgSLOPE, graphHUGE$omega))
        
        SNgLASSO <- c(SNgLASSO, SN(omegaHATgLASSO, graphHUGE$omega))
        SNgSLOPE <- c(SNgSLOPE, SN(omegaHATgLASSO, graphHUGE$omega))
        
        if(verbose)
            setTxtProgressBar(progressBar, i)
    }
    if(verbose)
        close(progressBar)
    
    return(list(gLASSO = list(FDR = sum(FDPgLASSO)/simulationsNumber,
                              Sensitivity = sum(SNgLASSO)/simulationsNumber),
                gSLOPE = list(FDR = sum(FDPgSLOPE)/simulationsNumber,
                              Sensitivity = sum(SNgSLOPE)/simulationsNumber)))
}