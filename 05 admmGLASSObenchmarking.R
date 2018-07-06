# Written by Micha³ Makowski

# install.packages("MASS")
# install.packages("glasso")
# install.packages("huge")

require(MASS)
require(glasso)
require(huge)
require(tictoc)

source("04 admmGLASSO.R")
source("01 auxilaryFunctions.R")



set.seed(100)

data  <- scale(dataSimulator(n = 100, SNR = 1, K = 5, numb.vars = 30, max.dim = 3)$X)
data = scale(huge.generator(n=400, d=100, graph="hub")$data) # Generate data with hub structures
data = scale(huge.generator(n=200, d=200, graph="cluster")$data) # Generate data with hub structures


tic("total")
for(mu in seq(from = .4, to = 1.5, length.out = 25))
{
    n     <- nrow(data)
    p     <- ncol(data)
    m     <- p*(p-1)/2
    alpha <- 0.05
    mBanerjee <- p^2   
    banerjeeLassoLambda <- qt(1-alpha/2/mBanerjee, df = n-2)/sqrt(n-2+qt(1-alpha/2/mBanerjee, df = n-2)^2)

    tic(mu)
    X <- glassoADMM(cov(data), mu = mu, lambda = banerjeeLassoLambda, penalizeDiagonal =  T)
    toc()
}
toc()

graphHUGE <- huge.generator(n = 100, d = 100, graph = "cluster")
precisionADMM <- $precisionMatrix
precisionGL   <- glasso(graphHUGE$sigmahat, rho = .034, penalize.diagonal = T)$wi

set.seed(100)
L = huge.generator(n=100,d=100,graph="hub") # Generate data with hub structures

X = L$data 
# X.pow = X^3/sqrt(15) 
X.npn = huge.npn(X.pow)

out = huge(X,nlambda=30., method = "glasso")
# out.mb = huge(X.pow,nlambda=30, method = "glasso") # Estimate the solution path
out.npn = huge(X.npn,nlambda=30, method = "glasso")

# huge.roc(out$path,L$theta)
# huge.roc(out$path,L$theta)
# huge.roc(out.mb$path,L$theta) # Plot the ROC curve
# huge.roc(out.npn$path,L$theta)

stars = huge.select(out,criterion="stars", stars.thresh=0.05) # Select the graph using StARS
# mb.stars = huge.select(out.mb,criterion="stars", stars.thresh=0.05) 
npn.stars = huge.select(out.npn,criterion="stars",stars.thresh=0.05)

ric = huge.select(out) # Select the graph using RIC
# mb.ric = huge.select(out.mb)
npn.ric = huge.select(out.npn)

precisionADMM <- glassoADMM(L$sigmahat, mu = .5, lambda = stars$opt.lambda, truncate = F)

plotCovariance(precisionADMM$precisionMatrix)
plotCovariance(as.matrix(stars$opt.icov))
plotCovariance(L$omega)

plotCovariance(precisionADMM$precisionMatrix - as.matrix(stars$opt.icov))


# plot(stars)
# plot(ric)
# plot(mb.ric)
# plot(npn.ric)
