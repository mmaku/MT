# Written by Micha³ Makowski

setwd("C:/Users/Martyna Strauchmann/Dropbox/04 Praca magisterska/gSLOPEtesting")

library(doSNOW)
library(foreach)
library(dplyr)
source("fastSimulationsFunctions.R")

# Cluster structure 
clusters <- createSimulationMatrix(nVec = c(50, 100, 200, 400),
                                   pVec = 100, 
                                   graphTypeVec = c("cluster"),
                                   alphaVec = c(0.05, 0.1, 0.2), 
                                   penalizeDiagonalVec = FALSE, 
                                   iterationsVec = 3000)

graphStructure_small_25 <- list(v = 0.3,
                                u = 0.7,
                                g = 10,
                                prob = 0.25)
graphStructure_small_50 <- list(v = 0.3,
                                u = 0.7,
                                g = 10,
                                prob = 0.5)
graphStructure_small_75 <- list(v = 0.3,
                                u = 0.7,
                                g = 10,
                                prob = 0.75)
graphStructure_small_100 <- list(v = 0.3,
                                 u = 0.7,
                                 g = 10,
                                 prob = 1)
graphStructure_large_25 <- list(v = 0.7,
                                u = 0.3,
                                g = 10,
                                prob = 0.25)
graphStructure_large_50 <- list(v = 0.7,
                                u = 0.3,
                                g = 10,
                                prob = 0.5)
graphStructure_large_75 <- list(v = 0.7,
                                u = 0.3,
                                g = 10,
                                prob = 0.75)
graphStructure_large_100 <- list(v = 0.7,
                                 u = 0.3,
                                 g = 10,
                                 prob = 1)

cl<-makeCluster(4) #your number of CPU cores
registerDoSNOW(cl)

doparList <-list(list(clusters, graphStructure_small_25), 
                 list(clusters, graphStructure_small_50), 
                 list(clusters, graphStructure_small_75), 
                 list(clusters, graphStructure_small_100), 
                 list(clusters, graphStructure_large_25), 
                 list(clusters, graphStructure_large_50), 
                 list(clusters, graphStructure_large_75), 
                 list(clusters, graphStructure_large_100))
results <- list() 
results <- foreach(i = doparList) %dopar% {
    source("fastSimulationsFunctions.R")
    simulations(i[[1]],
                saveAll = FALSE,
                graphParameters = i[[2]],
                fileName = paste0("_Prob_",i[[2]]$v,"_",i[[2]]$prob))
}

stopCluster(cl)

for(r in 1:length(results))
{
    graphParameters <- doparList[[r]][[2]]

    names(graphParameters) <- do.call(paste0, (expand.grid("graph.", names(graphParameters))))

    graphDF = data.frame(graphParameters)
    
    results[[r]] %>%
        merge(graphDF) -> results[[r]]
}


finalResults <- bind_rows(results)
filenameAll <- paste0("AllOneSimulations_Prob@", format(Sys.time(), '%y_%m_%d_%H_%M'), "#",nrow(finalResults))

save(finalResults, graphParameters,
     file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))


