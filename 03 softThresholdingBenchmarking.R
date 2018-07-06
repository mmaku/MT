# Written by Micha³ Makowski

# Comprasion of two diffrent implementation of soft thresholding function
# first using indexing, second using pmax

source("softThresholding.R")

# install.packages("tictoc")
require(tictoc)

tic("total")
tic("data generation")

set.seed(10)
testMatrices <- list()
resultsIndexing <- list()
resultsPmax <- list()

N <- 10000


for(i in 1:N)
    testMatrices[[i]] <- matrix(rnorm(400), 20, 20)

testThreshold <- rnorm(N, sd = .9)

toc()
tic("indexing")

for(i in 1:N)
    resultsIndexing[[i]] <- softThresholding(testMatrices[[i]], testThreshold[i])

toc()
tic("pmax")

for(i in 1:N)
    resultsPmax[[i]] <- softThresholding2(testMatrices[[i]], testThreshold[i])

toc()
toc()
