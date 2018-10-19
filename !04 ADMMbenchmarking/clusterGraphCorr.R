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

graphStructure_02_08 <- list(v = 0.2,
                             u = 0.8,
                             g = 10,
                             prob = 1)
graphStructure_03_07 <- list(v = 0.3,
                             u = 0.7,
                             g = 10,
                             prob = 1)
graphStructure_04_06 <- list(v = 0.4,
                             u = 0.6,
                             g = 10,
                             prob = 1)
graphStructure_06_04 <- list(v = 0.6,
                             u = 0.4,
                             g = 10,
                             prob = 1)
graphStructure_07_03 <- list(v = 0.7,
                             u = 0.3,
                             g = 10,
                             prob = 1)
graphStructure_08_02 <- list(v = 0.8,
                             u = 0.2,
                             g = 10,
                             prob = 1)


cl<-makeCluster(6) #your number of CPU cores
registerDoSNOW(cl)

doparList <-list(list(clusters, graphStructure_02_08), 
                 list(clusters, graphStructure_03_07), 
                 list(clusters, graphStructure_04_06), 
                 list(clusters, graphStructure_06_04), 
                 list(clusters, graphStructure_07_03), 
                 list(clusters, graphStructure_08_02))
results <- list() 
results <- foreach(i = doparList) %dopar% {
    source("fastSimulationsFunctions.R")
    simulations(i[[1]],
                saveAll = TRUE,
                testLocalFDR = FALSE,
                graphParameters = i[[2]],
                fileName = paste0("_Corr_",i[[2]]$v,"_",i[[2]]$u))
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
filenameAll <- paste0("AllOneSimulations_Corr@", format(Sys.time(), '%y_%m_%d_%H_%M'), "#",nrow(finalResults))

save(finalResults, graphParameters,
     file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))


