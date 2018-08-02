# Written by Micha³ Makowski

library(doSNOW)
library(foreach)
source("12 simulationsFunctions.R")

# iterationsMatrix <- createSimulationMatrix(nVec = 150, 
#                                            pVec = 200, 
#                                            graphTypeVec = c("hub", "cluster"),
#                                            alphaVec = c(0.1, 0.05, 0.01), 
#                                            penalizeDiagonalVec = FALSE, 
#                                            partialVec = FALSE,
#                                            iterationsVec = c(10, 20, 41, 83, 166, 312, 625, 1250, 2, 5000))

# iterationTest <- simulations(iterationsMatrix,
#                              saveAll = TRUE)


clusters <- createSimulationMatrix(nVec = c(100, 150, 200, 250, 300), 
                                            pVec = 200, 
                                            graphTypeVec = c("cluster"),
                                            alphaVec = c(0.05, 0.05*199/2/200), 
                                            penalizeDiagonalVec = FALSE, 
                                            partialVec = TRUE,
                                            iterationsVec = 2500)

# clustersTest <- simulations(clusters,
#                             saveAll = TRUE)

methods <- createSimulationMatrix(nVec = 200, 
                                  pVec = 200, 
                                  graphTypeVec = c("cluster"),
                                  alphaVec = c(0.05), 
                                  penalizeDiagonalVec = FALSE, 
                                  partialVec = TRUE,
                                  iterationsVec = 2500)

additional <- list(starsMB = list(method = "mb", criterion = "stars"),
                   ricMB = list(method = "mb"),
                   starsgLASSO = list(method = "glasso", criterion = "stars"),
                   ricgLASSO = list(method = "glasso"),
                   ebicgGLASSO = list(method = "glasso", criterion = "ebic"))

# methodsTest <- simulations(methods,
#                            saveAll = TRUE,
#                            additionalMethods = additional)

hubs <- createSimulationMatrix(nVec = c(50, 100, 150), 
                               pVec = 100, 
                               graphTypeVec = c("hub"),
                               alphaVec = c(0.05, 0.05*99/2/100), 
                               penalizeDiagonalVec = FALSE, 
                               partialVec = TRUE,
                               iterationsVec = 2500)

# hubsTest <- simulations(hubs,
#                         saveAll = TRUE)






cl<-makeCluster(3) #change the 2 to your number of CPU cores
registerDoSNOW(cl)

doparList <-list(list(clusters, NULL), list(methods, additional), list(hubs, NULL))
results <- list() 
results <- foreach(i = doparList) %dopar% {
    source("12 simulationsFunctions.R")
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
#                                      iterationsVec = 2500)
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