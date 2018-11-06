# Written by Micha³ Makowski

# setwd("C:/Users/Martyna Strauchmann/Dropbox/04 Praca magisterska/gSLOPEtesting")
# Start 11:43 06.11.2018
source("graphROCfunctions.R")

library(doSNOW)
library(foreach)

cl<-makeCluster(3) #your number of CPU cores
registerDoSNOW(cl)

doparList <-list(list("cluster", list(   v =  0.3, u = 0.7, prob = 0.5, g = 10)),
                 list("cluster", list(   v =  0.3, u = 0.7, prob = 1, g = 10)),
                 list("cluster", list(   v =  0.7, u = 0.3, prob = 0.5, g = 10)),
                 list("cluster", list(   v =  0.7, u = 0.3, prob = 1, g = 10)),
                 list("cluster", list(   v = -1,   u = 1,   prob = 0.5, g = 10)),
                 list("cluster", list(   v = -1,   u = 1,   prob = 1, g = 10)),
                 list("hub", list(       v =  0.3, u = 0.7, g = 10)),
                 list("hub", list(       v =  0.7, u = 0.3, g = 10)),
                 list("hub", list(       v = -1,   u = 1,   g = 10)),
                 list("scale-free", list(v =  0.3, u = 0.7, g = 10)),
                 list("scale-free", list(v =  0.7, u = 0.3, g = 10)),
                 list("scale-free", list(v = -1,   u = 1,   g = 10)))


resultsList <- list()
resultsList <- foreach(i = doparList) %dopar% {
    source("graphROCfunctions.R")
    list(graphType = i[[1]],
         graphParameters = i[[2]],
         results = graphROC(n = 100, 
                            p = 100, 
                            graphType = i[[1]],
                            graphParameters = i[[2]], 
                            alpha = 0.05, 
                            alphaMultiSeq = c(seq(from = 1e-4, 
                                                  to = 0.25, 
                                                  length.out = 15),
                                              seq(from = 0.25, 
                                                  to = 0.75, 
                                                  length.out = 16)[-1],
                                              seq(from = 0.75, 
                                                  to = 1.75, 
                                                  length.out = 16)[-1]),
                            iterations = 3000,
                            verbose = FALSE)) 
}

save(resultsList,
     file = paste0("./!02 Data/01 Binded/ROCfinal.RData"))
