# Written by Micha³ Makowski

# install.packages("glasso")
# install.packages("huge")

require(glasso)
require(huge)

source("11 measures.R")

createSimulationMatrix <- function(nVec = 150, 
                                   pVec = 200, 
                                   graphTypeVec = "cluster",
                                   alphaVec = 0.1, 
                                   penalizeDiagonalVec = FALSE, 
                                   partialVec = TRUE,
                                   iterationsVec = 1000)
{
    output <- expand.grid(nVec, pVec, graphTypeVec, alphaVec, penalizeDiagonalVec, partialVec, iterationsVec,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    colnames(output) <- c("n", "p", "graphType", "alpha", "penalizeDiagonal", "partial", "iterations")
    
    return(output)
}

simulations <- function(simulationMatrix, 
                        additionalMethods = NULL, 
                        verbose = TRUE,
                        saveEach = FALSE,
                        saveAll = FALSE)
{
    specificDoCall <- function(x) 
        doCall("measures", additionalMethods = additionalMethods, verbose = FALSE, args = x)
    
    if(verbose) 
    {
        cat("Starting simulations\nnumber of setups = ", NROW(simulationMatrix), 
            "\nnumber of simulations = ", sum(simulationMatrix$iterations), ".\n")
        
        progressBar <- txtProgressBar(min = 0, max = sum(simulationMatrix$iterations), style = 3)
        setTxtProgressBar(progressBar, 0)
    }
    if(saveAll)
    {
        filenameAll <- paste0("Simulation@", format(Sys.time(), '%y_%m_%d@%H_%M'), "#",
                              NROW(simulationMatrix)*(3+NROW(additionalMethods)))
    }
    
    output <- list()
    
    for(r in seq_len(NROW(simulationMatrix))) 
    {
        simResults <- specificDoCall(simulationMatrix[r,])
        
        output[[r]] <- cbind(simResults, simulationMatrix[r,], row.names = NULL)
        
        if(saveEach)
        {
            filename <- paste0("Simulation_", 
                               # format(Sys.time(), '%y_%m_%d_%H_%M'), 
                               paste(simulationMatrix[r,], collapse = '_'))
            
            setup <- simulationMatrix[r,]
            
            save(simResults, setup, additionalMethods, file = paste0("./!02 Data/", filename, ".RData"))
        }
        if(saveAll)
        {
            tempOutput <- do.call("rbind", output)
            save(tempOutput, additionalMethods, 
                 file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))
        }
        
        if(verbose)
            setTxtProgressBar(progressBar, sum(simulationMatrix$iterations[1:r]))
    }
    
    if(verbose)
        close(progressBar)
    
    output <- do.call("rbind", output)
    
    if(saveAll)
    {
        save(output, additionalMethods, 
             file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))
    }
    
    return(output)
}
