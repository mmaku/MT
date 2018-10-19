# Written by Micha³ Makowski

library(doSNOW)
library(foreach)
source("simulationsFunctions.R")

# Cluster structure 
clusterBHtest <- createSimulationMatrix(nVec = 75, 
                                        pVec = 100, 
                                        graphTypeVec = "cluster",
                                        alphaVec = c(0.05, 0.05*99/2/100), 
                                        penalizeDiagonalVec = FALSE, 
                                        iterationsVec = 3000)

cluster <- createSimulationMatrix(nVec = c(50, 100, 150),
                                  pVec = 100,
                                  graphTypeVec = "cluster",
                                  alphaVec = c(0.1, 0.05, 0.01),
                                  penalizeDiagonalVec = FALSE,
                                  iterationsVec = 3000)


# Scall-free structure14
scale3sBHtest <- createSimulationMatrix(nVec = 50,
                                 pVec = 50, 
                                 graphTypeVec = "hub",
                                 alphaVec = c(0.05, 0.05*49/2/50), 
                                 penalizeDiagonalVec = FALSE, 
                                 iterationsVec = 2000)

scale3 <- createSimulationMatrix(nVec = c(50, 100, 150),
                                 pVec = 100,
                                 graphTypeVec = "hub",
                                 alphaVec = c(0.1, 0.05, 0.01),
                                 penalizeDiagonalVec = FALSE,
                                 iterationsVec = 2000)

# Scall-free structure14
scale3sBHtest <- createSimulationMatrix(nVec = 50,
                                        pVec = 50, 
                                        graphTypeVec = "scale-free",
                                        alphaVec = c(0.05, 0.05*49/2/50), 
                                        penalizeDiagonalVec = FALSE, 
                                        iterationsVec = 3000)

scale3 <- createSimulationMatrix(nVec = c(50, 100, 150),
                                 pVec = 100,
                                 graphTypeVec = "scale-free",
                                 alphaVec = c(0.1, 0.05, 0.01),
                                 penalizeDiagonalVec = FALSE,
                                 iterationsVec = 3000)
# N increment test
Nincrement <- createSimulationMatrix(nVec = c(20,80,320,1280,5120),
                                     pVec = 80,
                                     graphTypeVec = "cluster",
                                     alphaVec = 0.05,
                                     penalizeDiagonalVec = FALSE,
                                     iterationsVec = 3000)


cl<-makeCluster(7) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

doparList <-list(list(clusters, NULL), 
                 list(scale3s, NULL),
                 list(hubs, NULL),
                 list(cluster, NULL),
                 list(scale3, NULL),
                 list(hub, NULL),
                 list(Nincrement, NULL))
results <- list() 
results <- foreach(i = doparList) %dopar% {
    source("simulationsFunctions.R")
    simulations(i[[1]],
                saveAll = TRUE,
                additionalMethods = i[[2]])
}

stopCluster(cl)

for(r in 1:length(results))
{
    filenameAll <- paste0("Simulation_", nrow(results[[r]]), "_",
                          format(Sys.time(), '%y_%m_%d_%H_%M'))
    
    output <- results[[r]]
    additionalMethods <- doparList[[r]][[2]]
    
    save(output, additionalMethods, 
         file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))
    
}




# 
# rocCluster <- createSimulationMatrix(nVec = 100, 
#                                      pVec = 100, 
#                                      graphTypeVec = c("cluster"),
#                                      alphaVec = seq(from = 0.01, to = 1, by = 0.01), 
#                                      penalizeDiagonalVec = FALSE, 
#                                      partialVec = TRUE,
#                                      iterationsVec = 3000)
# 
# rocClusterTest <- simulations(rocCluster,
#                               saveAll = FALSE)
# lambdaBanerjee(1000, 1000, alpha = 100)


# plot(1-rocClusterTest$SP[rocClusterTest$procedure == "gLASSO"], 
#      rocClusterTest$SN[rocClusterTest$procedure == "gLASSO"], "l")
# plot(1-rocClusterTest$SP[rocClusterTest$procedure == "holmgSLOPE"], 
#      rocClusterTest$SN[rocClusterTest$procedure == "holmgSLOPE"], "l")
# plot(1-rocClusterTest$SP[rocClusterTest$procedure == "BHgSLOPE"], 
#      rocClusterTest$SN[rocClusterTest$procedure == "BHgSLOPE"], "l")

# measures(n = 50, p = 100, graphType = "hub", alpha = 0.1, iterations = 100)


# nVec = 150, 
# pVec = 200, 
# graphTypeVec = "cluster",
# alphaVec = 0.05, 
# penalizeDiagonalVec = FALSE, 
# partialVec = FALSE,
# iterationsVec = 1000