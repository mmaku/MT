# Written by Micha³ Makowski

# setwd("C:/Users/Martyna Strauchmann/Dropbox/04 Praca magisterska/gSLOPEtesting")
# Start 11:43 06.11.2018
source("graphROCfunctions.R")

resultsList <- list()
iter <- 1

for(gt in c("cluster")) 
{
    for (gp in list(list(v = 0.3, u = 0.7, prob = 0.5, g = 10),
                    list(v = 0.3, u = 0.7, prob = 1, g = 10),
                    list(v = 0.7, u = 0.3, prob = 0.5, g = 10),
                    list(v = 0.7, u = 0.3, prob = 1, g = 10),
                    list(v = -1, u = 1, prob = 0.5, g = 10),
                    list(v = -1, u = 1, prob = 1, g = 10)))
    {
        if(!(gt %in% c("hub", "scale-free") & gp$prob == 1))
        {
            resultsList[[iter]] <- list(graphType = gt,
                                        graphParameters = gp,
                                        results = graphROC(n = 100, 
                                                           p = 100, 
                                                           graphType = gt,
                                                           graphParameters = gp, 
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
                                                           verbose = TRUE)) 
            iter <- iter + 1   
        }
    }
}


save(resultsList,
     file = paste0("./!02 Data/01 Binded/ROCfinal.RData"))

