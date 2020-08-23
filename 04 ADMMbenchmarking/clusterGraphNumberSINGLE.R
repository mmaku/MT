# Written by Micha³ Makowski

# setwd("C:/Users/Martyna Strauchmann/Dropbox/04 Praca magisterska/gSLOPEtesting")

library(dplyr)
source("fastSimulationsFunctions.R")

# Cluster structure 
clusters <- createSimulationMatrix(nVec = c(50, 100, 200, 400), 
                                   pVec = 100, 
                                   graphTypeVec = c("cluster", "hub"),
                                   alphaVec = c(0.05, 0.2), 
                                   penalizeDiagonalVec = FALSE, 
                                   iterationsVec = 3000)

graphStructure_equal_20 <- list(v = -1,
                                u = 1,
                                g = 20,
                                prob = 0.5)

results <- simulations(clusters,
                       saveAll = TRUE,
                       testLocalFDR = TRUE,
                       graphParameters = graphStructure_equal_20,
                       fileName = paste0("_HC_scaled_FDR_Num_BIS_",graphStructure_equal_20$v,"_",graphStructure_equal_20$g))

graphParameters <- graphStructure_equal_20

names(graphParameters) <- do.call(paste0, (expand.grid("graph.", names(graphParameters))))

graphDF = data.frame(graphParameters)

results %>%
    merge(graphDF) -> results

filenameAll <- paste0("AllOneSimulations_Num_BIS@", format(Sys.time(), '%y_%m_%d_%H_%M'), "#", nrow(results))

save(results, graphParameters,
     file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))


