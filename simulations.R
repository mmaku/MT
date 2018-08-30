# Written by Micha³ Makowski

library(doSNOW)
library(foreach)
source("simulationsFunctions.R")

# Cluster structure 
clusters <- createSimulationMatrix(nVec = c(50, 100, 150), 
                                   pVec = 100, 
                                   graphTypeVec = c("cluster"),
                                   alphaVec = c(0.05, 0.05*99/2/100), 
                                   penalizeDiagonalVec = FALSE, 
                                   iterationsVec = 3000)

# Scall-free structure
scale3 <- createSimulationMatrix(nVec = c(25, 50, 75),
                                 pVec = 50, 
                                 graphTypeVec = c("scale-free"),
                                 alphaVec = c(0.05, 0.05*49/2/50), 
                                 penalizeDiagonalVec = FALSE, 
                                 iterationsVec = 3000)

# Hub structure
hubs <- createSimulationMatrix(nVec = c(100, 200, 300),
                               pVec = 200, 
                               graphTypeVec = c("cluster"),
                               alphaVec = c(0.05, 0.05*199/2/200), 
                               penalizeDiagonalVec = FALSE, 
                               iterationsVec = 3000)

main <- createSimulationMatrix(nVec = 75,
                               pVec = 100,
                               graphTypeVec = c("hub", "cluster"),
                               alphaVec = c(0.1, 0.05, 0.01),
                               penalizeDiagonalVec = FALSE,
                               iterationsVec = 3000)


# Methods comparison
# methods <- createSimulationMatrix(nVec = 200, 
#                                   pVec = 200, 
#                                   graphTypeVec = c("cluster"),
#                                   alphaVec = c(0.05), 
#                                   penalizeDiagonalVec = FALSE, 
#                                   iterationsVec = 3000)
# 
# additional <- list(starsMB = list(method = "mb", criterion = "stars"),
#                    ricMB = list(method = "mb"),
#                    starsgLASSO = list(method = "glasso", criterion = "stars"),
#                    ricgLASSO = list(method = "glasso"),
#                    ebicgGLASSO = list(method = "glasso", criterion = "ebic"))

# methodsTest <- simulations(methods,
#                            saveAll = TRUE,
#                            additionalMethods = additional)


cl<-makeCluster(4) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

doparList <-list(list(clusters, NULL), 
                 list(scale3, NULL),
                 list(hubs, NULL),
                 # list(methods, additional),
                 list(main, NULL))
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