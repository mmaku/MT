# Written by Micha³ Makowski

# setwd("C:/Users/Martyna Strauchmann/Dropbox/04 Praca magisterska/gSLOPEtesting")

source("graphROCfunctions.R")

resultsList <- list()
iter <- 1

for(gt in c("cluster", "hub", "scale-free")) 
{
    for (gp in list(list(v = 0.3, u = 0.7, prob = 0.5, g = 10), 
                    list(v = 0.7, u = 0.3, prob = 0.5, g = 10),
                    list(v = -1, u = 1, prob = 0.5, g = 10))) 
    {
        resultsList[[iter]] <- list(graphType = gt,
                                    graphParameters = gp,
                                    results = graphROC(n = 100, 
                                                       p = 100, 
                                                       graphType = gt,
                                                       graphParameters = gp, 
                                                       alpha = 0.05, 
                                                       alphaMultiSeq = c(seq(from = 1e-3, 
                                                                             to = 1.5, 
                                                                             length.out = 30)),
                                                       iterations = 1000,
                                                       verbose = FALSE)) 
        iter <- iter + 1
    }
}


save(resultsList,
     file = paste0("./!02 Data/01 Binded/ROC3.RData"))


# cl<-makeCluster(2) #your number of CPU cores
# registerDoSNOW(cl)
# 
# doparList <-list(list(clusters, graphStructure_small), 
#                  list(hubs, graphStructure_small), 
#                  list(scale3, graphStructure_small), 
#                  list(clusters, graphStructure_large), 
#                  list(hubs, graphStructure_large), 
#                  list(scale3, graphStructure_large))
# 
# results <- list() 
# results <- foreach(i = doparList) %dopar% {
#     source("fastSimulationsFunctions.R")
#     graphROC(i[[1]],
#                 saveAll = TRUE,
#                 testLocalFDR = TRUE,
#                 graphParameters = i[[2]],
#                 fileName = paste0("_Cl_both_FDR_Prob_",i[[2]]$v,"_",i[[2]]$prob))
# }
# 
# stopCluster(cl)
# 
# for(r in 1:length(results))
# {
#     graphParameters <- doparList[[r]][[2]]
# 
#     names(graphParameters) <- do.call(paste0, (expand.grid("graph.", names(graphParameters))))
# 
#     graphDF = data.frame(graphParameters)
#     
#     results[[r]] %>%
#         merge(graphDF) -> results[[r]]
# }
# 
# 
# finalResults <- bind_rows(results)
# filenameAll <- paste0("AllOneSimulations_Cl_both_FDR_Prob@", format(Sys.time(), '%y_%m_%d_%H_%M'), "#",nrow(finalResults))
# 
# save(finalResults, graphParameters,
#      file = paste0("./!02 Data/01 Binded/", filenameAll, ".RData"))
# 
# 
