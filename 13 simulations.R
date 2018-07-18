# Written by Micha³ Makowski

source("12 simulationsFunctions.R")

# iterationsMatrix <- createSimulationMatrix(nVec = 150, 
#                                            pVec = 200, 
#                                            graphTypeVec = c("hub", "cluster"),
#                                            alphaVec = c(0.1, 0.05, 0.01), 
#                                            penalizeDiagonalVec = FALSE, 
#                                            partialVec = FALSE,
#                                            iterationsVec = c(10, 20, 41, 83, 166, 312, 625, 1250, 2500, 5000))

# iterationTest <- simulations(iterationsMatrix,
#                              saveAll = TRUE)


sanityCheckMatrix <- createSimulationMatrix(nVec = c(100, 150, 200, 250, 300), 
                                            pVec = 200, 
                                            graphTypeVec = c("cluster"),
                                            alphaVec = c(0.05, 0.05*199/2/200), 
                                            penalizeDiagonalVec = FALSE, 
                                            partialVec = TRUE,
                                            iterationsVec = 250)

sanityCheck <- simulations(sanityCheckMatrix,
                           saveAll = TRUE)


firstTestMatrix <- createSimulationMatrix(nVec = c(100, 150, 200, 250, 300), 
                                          pVec = 200, 
                                          graphTypeVec = c("hub", "cluster"),
                                          alphaVec = c(0.1, 0.05, 0.01), 
                                          penalizeDiagonalVec = FALSE, 
                                          partialVec = TRUE,
                                          iterationsVec = 250)

firstTest <- simulations(firstTestMatrix,
                         saveAll = TRUE)

# additional <- list(starsMB = list(method = "mb", criterion = "stars"),
#                    ricMB = list(method = "mb"),
#                    starsgLASSO = list(method = "glasso", criterion = "stars"),
#                    ricgLASSO = list(method = "glasso"),
#                    ebicgGLASSO = list(method = "glasso", criterion = "ebic"))

# nVec = 150, 
# pVec = 200, 
# graphTypeVec = "cluster",
# alphaVec = 0.05, 
# penalizeDiagonalVec = FALSE, 
# partialVec = FALSE,
# iterationsVec = 1000