# Written by Micha³ Makowski

# install.packages("glasso")
# install.packages("huge")
# install.packages("microbenmark")

# devtools::install_github("olafmersmann/microbenchmarkCore")
# devtools::install_github("olafmersmann/microbenchmark")

require(glasso)
require(huge)
require(microbenchmark)
require(xtable)

source("lambdas.R")
source("auxilaryFunctions.R")
source("admmGLASSO.R")


set.seed(100)

hubData <- scale(huge.generator(n=400, d=100, graph="hub")$data) # Generate data with hub structures
clusterData <- scale(huge.generator(n=200, d=200, graph="cluster")$data) # Generate data with hub structures

hubLambda <- lambdaSelector(hubData, method="banerjee")
clusterLambda <- lambdaSelector(clusterData, method="banerjee")

hubCov <- cov(hubData)
clusterCov <- cov(clusterData)

# results <- matrix(0, nrow=length(muSeries))

# Hub
a <- microbenchmark(glassoADMM_hub=glassoADMM(hubCov, lambda=hubLambda, penalizeDiagonal=FALSE, verbose=FALSE),
                           glasso_hub=glasso(hubCov, rho=hubLambda, penalize.diagonal=FALSE))

hubPrecADMM <- glassoADMM(hubCov, lambda=hubLambda, penalizeDiagonal=FALSE, verbose=FALSE)$precisionMatrix
hubPrecGL <- glasso(hubCov, rho=hubLambda, penalize.diagonal=FALSE)$wi

hubObjectiveADMM <- -log(det(hubPrecADMM)) + sum(diag(hubPrecADMM%*%hubCov)) + hubLambda*sum(abs(hubPrecADMM))
hubObjectiveGL <- -log(det(hubPrecGL)) + sum(diag(hubPrecGL%*%hubCov)) + hubLambda*sum(abs(hubPrecGL))

# Cluster
clusterBench <- microbenchmark(glassoADMM_cluster=glassoADMM(clusterCov, lambda=clusterLambda, penalizeDiagonal=FALSE, verbose=FALSE),
                               glasso_cluster=glasso(clusterCov, rho=clusterLambda, penalize.diagonal=FALSE))

clusterPrecADMM <- glassoADMM(clusterCov, lambda=clusterLambda, penalizeDiagonal=FALSE, verbose=FALSE)$precisionMatrix
clusterPrecGL <- glasso(clusterCov, rho=clusterLambda, penalize.diagonal=FALSE)$wi

clusterObjectiveADMM <- -log(det(clusterPrecADMM)) + sum(diag(clusterPrecADMM%*%clusterCov)) + clusterLambda*sum(abs(clusterPrecADMM))
clusterObjectiveGL <- -log(det(clusterPrecGL)) + sum(diag(clusterPrecGL%*%clusterCov)) + clusterLambda*sum(abs(clusterPrecGL))

# Results tables

italic <- function(x)
{
    paste0('{\\textit{', x, '}}')
}

hubTiming <- summary(hubBench)[,c(-1,-3,-6,-9)]
hubTiming <-  round(hubTiming, 4)
rownames(hubTiming) <- c("ADMM", "graphical Lasso")
colnames(hubTiming) <- c("Minimum", "Mean", "Median", "Maximum", "# of evaluations")

clusterTiming <- summary(clusterBench)[,c(-1,-3,-6,-9)]
clusterTiming <-  round(clusterTiming, 4)
rownames(clusterTiming) <- c("ADMM", "graphical Lasso")
colnames(clusterTiming) <- c("Minimum", "Mean", "Median", "Maximum", "# of evaluations")

objectiveFunction <- data.frame(list(hub=c(hubObjectiveADMM, hubObjectiveGL), 
                                     cluster = c(clusterObjectiveADMM, clusterObjectiveGL)),
                                row.names=c("ADMM", "graphical Lasso"))
objectiveFunction <-  round(objectiveFunction, 4)


print(xtable(hubTiming, 
             caption = "Algorithms execution times comparison for the hub graph structure.",
             auto = TRUE),
      sanitize.rownames.function = italic,
      booktabs = TRUE)

print(xtable(clusterTiming, 
             caption = "Algorithms execution times comparison for the cluster graph structure.",
             auto = TRUE),
      booktabs = TRUE)


print(xtable(objectiveFunction, 
             caption = "Objective function value for optimal point found by each algoritmh in two different settings.",
             auto = TRUE),
      booktabs = TRUE)



