# Written by Micha³ Makowski

# install.packages("glasso")
# install.packages("huge")

require(glasso)
require(huge)
require(tictoc)

source("fastMeasures.R")

createSimulationMatrix <- function(nVec = 150, 
                                   pVec = 200, 
                                   graphTypeVec = "cluster",
                                   alphaVec = 0.05,
                                   penalizeDiagonalVec = FALSE, 
                                   iterationsVec = 1000)
{
    output <- expand.grid(nVec, pVec, graphTypeVec, alphaVec, penalizeDiagonalVec, iterationsVec,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    colnames(output) <- c("n", "p", "graphType", "alpha", "penalizeDiagonal", "iterations")
    
    return(output)
}

simulations <- function(simulationMatrix, 
                        graphParameters = NULL, 
                        saveAll = FALSE,
                        fileName = "",
                        testLocalFDR = TRUE)
{
    tic("All")
    ticFile <- paste0("./!02 Data/01 Binded/tic", fileName, ".txt")
    
    specificDoCall <- function(x) 
        doCall("measures", 
               graphParameters = graphParameters, testLocalFDR = testLocalFDR, verbose = FALSE, 
               args = x)
    
    if(saveAll)
    {
        filenameAll <- paste0("AllSimulations", fileName, "@", format(Sys.time(), '%y_%m_%d@%H_%M'), "#",
                              NROW(simulationMatrix)*3)
    }
    
    output <- list()
    
    for(r in seq_len(NROW(simulationMatrix))) 
    {
        ticName <- paste(simulationMatrix[r,], collapse = '_')
        tic(ticName)
        
        simResults <- specificDoCall(simulationMatrix[r,])
        output[[r]] <- cbind(simResults, simulationMatrix[r,], row.names = NULL)
        
        if(saveAll)
        {
            tempOutput <- do.call("rbind", output)
            save(tempOutput, graphParameters, 
                 file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))
        }
        
        toc(quiet = TRUE, log = TRUE)
        write(tic.log()[[1]], file = ticFile, append = TRUE)
        tic.clearlog()
    }
    
    output <- do.call("rbind", output)
    
    if(saveAll)
    {
        save(output, graphParameters, 
             file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))
    }
    
    toc(quiet = TRUE, log = TRUE)
    write(tic.log()[[1]], file = ticFile, append = TRUE)
    tic.clearlog()
    
    return(output)
}
