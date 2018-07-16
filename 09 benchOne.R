# Written by Micha³ Makowski

# install.packages("MASS")
# install.packages("glasso")
# install.packages("huge")

require(MASS)
require(glasso)
require(huge)

source("01 auxilaryFunctions.R")
source("04 admmGLASSO.R")
source("07 admmGSLOPE.R")
source("08 benchOneFunctions")


set.seed(100)

data  <- scale(dataSimulator(n = 100, SNR = 1, K = 5, numb.vars = 30, max.dim = 3)$X)

str(data)

firstBench_0.1_F_T <- benchOne(data, alpha = .1)
firstBench_0.1_F_F <- benchOne(data, alpha = .1, truncate = F)
firstBench_0.1_T_T <- benchOne(data, alpha = .1, penalizeDiagonal = T)
firstBench_0.1_T_F <- benchOne(data, alpha = .1, truncate = F, penalizeDiagonal = T)

firstBench_0.05_F_T <- benchOne(data)
firstBench_0.05_F_F <- benchOne(data, truncate = F)
firstBench_0.05_T_T <- benchOne(data, penalizeDiagonal = T)
firstBench_0.05_T_F <- benchOne(data, truncate = F, penalizeDiagonal = T)

firstBench_0.01_F_T <- benchOne(data, alpha = .01)
firstBench_0.01_F_F <- benchOne(data, alpha = .01, truncate = F)
firstBench_0.01_T_T <- benchOne(data, alpha = .01, penalizeDiagonal = T)
firstBench_0.01_T_F <- benchOne(data, alpha = .01, truncate = F, penalizeDiagonal = T)

firstBench_0.05_F_T_ric <- benchOne(data, selectCriterion = "ric")
firstBench_0.05_F_T_ebic <- benchOne(data, selectCriterion = "ebic")

ggsave("firstBench_0.05_F_T.png", plotBenchOne(firstBench_0.05_F_T), device = "png", path = "./!01 Plots/01 first/01 criterion/", dpi = 200)
ggsave("firstBench_0.05_F_T_ric.png", plotBenchOne(firstBench_0.05_F_T_ric), device = "png", path = "./!01 Plots/01 first/01 criterion/", dpi = 200)
ggsave("firstBench_0.05_F_T_ebic.png", plotBenchOne(firstBench_0.05_F_T_ebic), device = "png", path = "./!01 Plots/01 first/01 criterion/", dpi = 200)


ggsave("firstBench_0.1_F_T.png", plotBenchOne(firstBench_0.1_F_T), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.1_T_F.png", plotBenchOne(firstBench_0.1_T_F), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.1_F_F.png", plotBenchOne(firstBench_0.1_F_F), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.1_T_T.png", plotBenchOne(firstBench_0.1_T_T), device = "png", path = "./!01 Plots/01 first/", dpi = 200)

ggsave("firstBench_0.05_F_T.png", plotBenchOne(firstBench_0.05_F_T), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.05_T_F.png", plotBenchOne(firstBench_0.05_T_F), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.05_F_F.png", plotBenchOne(firstBench_0.05_F_F), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.05_T_T.png", plotBenchOne(firstBench_0.05_T_T), device = "png", path = "./!01 Plots/01 first/", dpi = 200)

ggsave("firstBench_0.01_F_T.png", plotBenchOne(firstBench_0.01_F_T), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.01_T_F.png", plotBenchOne(firstBench_0.01_T_F), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.01_F_F.png", plotBenchOne(firstBench_0.01_F_F), device = "png", path = "./!01 Plots/01 first/", dpi = 200)
ggsave("firstBench_0.01_T_T.png", plotBenchOne(firstBench_0.01_T_T), device = "png", path = "./!01 Plots/01 first/", dpi = 200)

##################################################################################################

data = scale(huge.generator(n=50, d=75, graph="hub")$data) # Generate data with hub structures

secondHubBench_0.1_F_T <- benchOne(data, alpha = .1)
secondHubBench_0.1_F_F <- benchOne(data, alpha = .1, truncate = F)
secondHubBench_0.1_T_T <- benchOne(data, alpha = .1, penalizeDiagonal = T)
secondHubBench_0.1_T_F <- benchOne(data, alpha = .1, truncate = F, penalizeDiagonal = T)

secondHubBench_0.05_F_T <- benchOne(data)
secondHubBench_0.05_F_F <- benchOne(data, truncate = F)
secondHubBench_0.05_T_T <- benchOne(data, penalizeDiagonal = T)
secondHubBench_0.05_T_F <- benchOne(data, truncate = F, penalizeDiagonal = T)

secondHubBench_0.01_F_T <- benchOne(data, alpha = .01)
secondHubBench_0.01_F_F <- benchOne(data, alpha = .01, truncate = F)
secondHubBench_0.01_T_T <- benchOne(data, alpha = .01, penalizeDiagonal = T)
secondHubBench_0.01_T_F <- benchOne(data, alpha = .01, truncate = F, penalizeDiagonal = T)

ggsave("secondHubBench_0.1_F_T.png", plotBenchOne(secondHubBench_0.1_F_T), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.1_T_F.png", plotBenchOne(secondHubBench_0.1_T_F), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.1_F_F.png", plotBenchOne(secondHubBench_0.1_F_F), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.1_T_T.png", plotBenchOne(secondHubBench_0.1_T_T), device = "png", path = "./!01 Plots/02 second/", dpi = 200)

ggsave("secondHubBench_0.05_F_T.png", plotBenchOne(secondHubBench_0.05_F_T), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.05_T_F.png", plotBenchOne(secondHubBench_0.05_T_F), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.05_F_F.png", plotBenchOne(secondHubBench_0.05_F_F), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.05_T_T.png", plotBenchOne(secondHubBench_0.05_T_T), device = "png", path = "./!01 Plots/02 second/", dpi = 200)

ggsave("secondHubBench_0.01_F_T.png", plotBenchOne(secondHubBench_0.01_F_T), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.01_T_F.png", plotBenchOne(secondHubBench_0.01_T_F), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.01_F_F.png", plotBenchOne(secondHubBench_0.01_F_F), device = "png", path = "./!01 Plots/02 second/", dpi = 200)
ggsave("secondHubBench_0.01_T_T.png", plotBenchOne(secondHubBench_0.01_T_T), device = "png", path = "./!01 Plots/02 second/", dpi = 200)

#######################################

data = scale(huge.generator(n=200, d=200, graph="cluster")$data) # Generate data with hub structures

thirdClusterBench_0.1_F_T <- benchOne(data, alpha = .1)
thirdClusterBench_0.1_F_F <- benchOne(data, alpha = .1, truncate = F)
thirdClusterBench_0.1_T_T <- benchOne(data, alpha = .1, penalizeDiagonal = T)
thirdClusterBench_0.1_T_F <- benchOne(data, alpha = .1, truncate = F, penalizeDiagonal = T)

thirdClusterBench_0.05_F_T <- benchOne(data)
thirdClusterBench_0.05_F_F <- benchOne(data, truncate = F)
thirdClusterBench_0.05_T_T <- benchOne(data, penalizeDiagonal = T)
thirdClusterBench_0.05_T_F <- benchOne(data, truncate = F, penalizeDiagonal = T)

thirdClusterBench_0.01_F_T <- benchOne(data, alpha = .01)
thirdClusterBench_0.01_F_F <- benchOne(data, alpha = .01, truncate = F)
thirdClusterBench_0.01_T_T <- benchOne(data, alpha = .01, penalizeDiagonal = T)
thirdClusterBench_0.01_T_F <- benchOne(data, alpha = .01, truncate = F, penalizeDiagonal = T)

ggsave("thirdClusterBench_0.1_F_T.png", plotBenchOne(thirdClusterBench_0.1_F_T), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.1_T_F.png", plotBenchOne(thirdClusterBench_0.1_T_F), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.1_F_F.png", plotBenchOne(thirdClusterBench_0.1_F_F), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.1_T_T.png", plotBenchOne(thirdClusterBench_0.1_T_T), device = "png", path = "./!01 Plots/03 third/", dpi = 200)

ggsave("thirdClusterBench_0.05_F_T.png", plotBenchOne(thirdClusterBench_0.05_F_T), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.05_T_F.png", plotBenchOne(thirdClusterBench_0.05_T_F), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.05_F_F.png", plotBenchOne(thirdClusterBench_0.05_F_F), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.05_T_T.png", plotBenchOne(thirdClusterBench_0.05_T_T), device = "png", path = "./!01 Plots/03 third/", dpi = 200)

ggsave("thirdClusterBench_0.01_F_T.png", plotBenchOne(thirdClusterBench_0.01_F_T), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.01_T_F.png", plotBenchOne(thirdClusterBench_0.01_T_F), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.01_F_F.png", plotBenchOne(thirdClusterBench_0.01_F_F), device = "png", path = "./!01 Plots/03 third/", dpi = 200)
ggsave("thirdClusterBench_0.01_T_T.png", plotBenchOne(thirdClusterBench_0.01_T_T), device = "png", path = "./!01 Plots/03 third/", dpi = 200)              
