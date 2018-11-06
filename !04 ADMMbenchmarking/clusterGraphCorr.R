# Written by Micha³ Makowski

setwd("C:/Users/Martyna Strauchmann/Dropbox/04 Praca magisterska/gSLOPEtesting")

library(doSNOW)
library(foreach)
library(dplyr)
source("fastSimulationsFunctions.R")

# Cluster structure 
clusters <- createSimulationMatrix(nVec = c(50, 100, 200, 400), 
                                   pVec = 100, 
                                   graphTypeVec = c("cluster", "hub", "scale-free"),
                                   alphaVec = c(0.05, 0.2),
                                   scaledVec = TRUE, 
                                   iterationsVec = 3000) 

graphStructure_07_03 <- list(v = 0.7,
                             u = 0.3,
                             g = 10,
                             prob = 0.5)
graphStructure_08_02 <- list(v = 0.8,
                             u = 0.2,
                             g = 10,
                             prob = 0.5)
graphStructure_m05_1 <- list(v = -0.5,
                             u = 1,
                             g = 10,
                             prob = 0.5)
graphStructure_m1_1 <- list(v = -1,
                             u = 1,
                             g = 10,
                             prob = 0.5)

cl<-makeCluster(2) #your number of CPU cores
registerDoSNOW(cl)

doparList <-list(list(clusters, graphStructure_07_03),
                 list(clusters, graphStructure_08_02),
                 list(clusters, graphStructure_m05_1),
                 list(clusters, graphStructure_m1_1))

results <- list() 
results <- foreach(i = doparList) %dopar% {
    source("fastSimulationsFunctions.R")
    simulations(i[[1]],
                saveAll = TRUE,
                testLocalFDR = TRUE,
                graphParameters = i[[2]],
                fileName = paste0("_HCSF_scaled_FDR_Corr_BIS_",i[[2]]$v,"_",i[[2]]$u))
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
filenameAll <- paste0("AllOneSimulations_HCSF_scaled_FDR_Corr_BIS@", format(Sys.time(), '%y_%m_%d_%H_%M'), "#",nrow(finalResults))

save(finalResults, graphParameters,
     file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))


