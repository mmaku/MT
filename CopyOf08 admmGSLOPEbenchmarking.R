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


set.seed(100)

data  <- scale(dataSimulator(n = 100, SNR = 1, K = 5, numb.vars = 30, max.dim = 3)$X)
sampleCovariance <- cov(data)

n     <- nrow(data)
p     <- ncol(data)
m     <- p*(p-1)/2
alpha <- 0.05
m_banerjee <- p^2
lambda <- qt(1-alpha/2/m_banerjee, df = n-2)/sqrt(n-2+qt(1-alpha/2/m_banerjee, df = n-2)^2)
precisionGLASSO <- glasso(sampleCovariance, rho = lambda)$wi
precisionGLASSOADMM <- glassoADMM(sampleCovariance, lambda = lambda, penalizeDiagonal = T, truncate = T)$precisionMatrix

k = 1:m
lambdaBis <- qt(1-alpha*k/2/m, df = n-2)/sqrt(n-2+qt(1-alpha*k/2/m, df = n-2)^2)
lambdaSeries <- c(rep(lambdaBis[1], p), rep(lambdaBis, each=2))
precisionGSLOPEADMM <- gslopeADMM(sampleCovariance, lambda = lambdaSeries, penalizeDiagonal = T, truncate = T)$precisionMatrix

alphaSmall <- 0.01
lambdaCis <- qt(1-alphaSmall*k/2/m, df = n-2)/sqrt(n-2+qt(1-alphaSmall*k/2/m, df = n-2)^2)
lambdaSeriesCis <- c(rep(lambdaCis[1], p), rep(lambdaCis, each=2))
precisionGSLOPEADMMsmallAlpha <- gslopeADMM(sampleCovariance, lambda = lambdaSeriesCis, penalizeDiagonal = F, truncate = T)$precisionMatrix

glassoHuge = huge(data, nlambda=30, method = "glasso")
astars = huge.select(glassoHuge, criterion="stars", stars.thresh=0.05) # Select the graph using StARS
aprecisionGLASSOADMMstars <- glassoADMM(sampleCovariance, mu = .5, lambda = astars$opt.lambda, penalizeDiagonal = F, truncate = T)$precisionMatrix

plotMatrix(aprecisionGLASSOADMMstars)
plotMatrix(as.matrix(astars$opt.icov))

plotMatrix(precisionGLASSO)
plotMatrix(precisionGLASSOADMM)
plotMatrix(precisionGSLOPEADMM)
plotMatrix(precisionGSLOPEADMMsmallAlpha)

##################################################################################################

L = huge.generator(n=200,d=50,graph="hub") # Generate data with hub structures

X = scale(L$data) 
X.pow = X^3/sqrt(15)
X.npn = huge.npn(X.pow)

out = huge(X,nlambda=30., method = "glasso")
out.mb = huge(X.pow,nlambda=30, method = "glasso") # Estimate the solution path
out.npn = huge(X.npn,nlambda=30, method = "glasso")

stars = huge.select(out,criterion="stars", stars.thresh=0.05) # Select the graph using StARS
mb.stars = huge.select(out.mb,criterion="stars", stars.thresh=0.05)
npn.stars = huge.select(out.npn,criterion="stars",stars.thresh=0.05)

n     <- nrow(L$data)
p     <- ncol(L$data)
m     <- p*(p-1)/2
alpha <- 0.05

k = 1:m
lambdaBis <- qt(1-alpha*k/2/m, df = n-2)/sqrt(n-2+qt(1-alpha*k/2/m, df = n-2)^2)
lambdaSeries <- c(rep(lambdaBis[1], p), rep(lambdaBis, each=2))
precisionGSLOPEADMM.hug <- gslopeADMM(cov(X), lambda = lambdaSeries)$precisionMatrix
precisionGSLOPEADMM.pow <- gslopeADMM(cov(X.pow), lambda = lambdaSeries)$precisionMatrix
precisionGSLOPEADMM.npn <- gslopeADMM(cov(X.npn), lambda = lambdaSeries)$precisionMatrix

alphaSmall <- 0.01
lambdaCis <- qt(1-alphaSmall*k/2/m, df = n-2)/sqrt(n-2+qt(1-alphaSmall*k/2/m, df = n-2)^2)
lambdaSeriesCis <- c(rep(lambdaCis[1], p), rep(lambdaCis, each=2))
precisionGSLOPEADMMsmallAlpha.hug <- gslopeADMM(cov(X), lambda = lambdaSeriesCis)$precisionMatrix
precisionGSLOPEADMMsmallAlpha.pow <- gslopeADMM(cov(X.pow), lambda = lambdaSeriesCis)$precisionMatrix
precisionGSLOPEADMMsmallAlpha.npn <- gslopeADMM(cov(X.npn), lambda = lambdaSeriesCis)$precisionMatrix

plotMatrix(precisionGSLOPEADMM.hug)
plotMatrix(precisionGSLOPEADMM.pow)
plotMatrix(precisionGSLOPEADMM.npn)
plotMatrix(precisionGSLOPEADMMsmallAlpha.hug)
plotMatrix(precisionGSLOPEADMMsmallAlpha.pow)
plotMatrix(precisionGSLOPEADMMsmallAlpha.npn)
plotMatrix(as.matrix(stars$opt.icov))
plotMatrix(as.matrix(mb.stars$opt.icov))
plotMatrix(as.matrix(npn.stars$opt.icov))
plotMatrix(solve(cov(X)))
plotMatrix(solve(cov(X.pow)))
plotMatrix(solve(cov(X.npn)))



#######################################

clL = huge.generator(n=200,d=50,graph="cluster") # Generate data with hub structures

clX = scale(clL$data) 
clX.pow = clX^3/sqrt(15)
clX.npn = huge.npn(clX.pow)

clout = huge(clX,nlambda=30., method = "glasso")
clout.mb = huge(clX.pow,nlambda=30, method = "glasso") # Estimate the solution path
clout.npn = huge(clX.npn,nlambda=30, method = "glasso")

clstars = huge.select(clout,criterion="stars", stars.thresh=0.05) # Select the graph using StARS
clmb.stars = huge.select(clout.mb,criterion="stars", stars.thresh=0.05)
clnpn.stars = huge.select(clout.npn,criterion="stars",stars.thresh=0.05)

cln     <- nrow(L$data)
clp     <- ncol(L$data)
clm     <- p*(p-1)/2
clalpha <- 0.05

clk = 1:clm
cllambdaBis <- qt(1-clalpha*clk/2/clm, df = cln-2)/sqrt(cln-2+qt(1-clalpha*clk/2/clm, df = cln-2)^2)
cllambdaSeries <- c(rep(cllambdaBis[1], clp), rep(cllambdaBis, each=2))
clprecisionGSLOPEADMM.hug <- gslopeADMM(cov(clX), lambda = cllambdaSeries)$precisionMatrix
clprecisionGSLOPEADMM.pow <- gslopeADMM(cov(clX.pow), lambda = cllambdaSeries)$precisionMatrix
clprecisionGSLOPEADMM.npn <- gslopeADMM(cov(clX.npn), lambda = cllambdaSeries)$precisionMatrix

alphaSmall <- 0.01
lambdaCis <- qt(1-alphaSmall*k/2/m, df = n-2)/sqrt(n-2+qt(1-alphaSmall*k/2/m, df = n-2)^2)
lambdaSeriesCis <- c(rep(lambdaCis[1], p), rep(lambdaCis, each=2))
clprecisionGSLOPEADMMsmallAlpha.hug <- gslopeADMM(cov(clX), lambda = lambdaSeriesCis)$precisionMatrix
clprecisionGSLOPEADMMsmallAlpha.pow <- gslopeADMM(cov(clX.pow), lambda = lambdaSeriesCis)$precisionMatrix
clprecisionGSLOPEADMMsmallAlpha.npn <- gslopeADMM(cov(clX.npn), lambda = lambdaSeriesCis)$precisionMatrix

plotMatrix(clprecisionGSLOPEADMM.hug)
plotMatrix(clprecisionGSLOPEADMM.pow)
plotMatrix(clprecisionGSLOPEADMM.npn)
plotMatrix(clprecisionGSLOPEADMMsmallAlpha.hug)
plotMatrix(clprecisionGSLOPEADMMsmallAlpha.pow)
plotMatrix(precisionGSLOPEADMMsmallAlpha.npn)
plotMatrix(as.matrix(clstars$opt.icov))
plotMatrix(as.matrix(clmb.stars$opt.icov))
plotMatrix(as.matrix(clnpn.stars$opt.icov))
plotMatrix(solve(cov(clX)))
plotMatrix(clL$omega)
plotMatrix(solve(cov(clX.pow)))
plotMatrix(solve(cov(clX.npn)))



#######################################


# data(stockdata)                                           # Load the stock data
# Y = log(stockdata$data[2:1258,]/stockdata$data[1:1257,])  # Preprocessing
# Y.npn = huge.npn(Y, npn.func="truncation")                # Nonparanormal
# out.npn = huge(Y.npn,method = "glasso", nlambda=40,lambda.min.ratio = 0.4)
# out = huge(Y,method = "glasso", nlambda=40,lambda.min.ratio = 0.4)
# 
# 
# n     <- nrow(Y)
# p     <- ncol(Y)
# m     <- p*(p-1)/2
# alpha <- 0.05
# 
# k = 1:m
# lambdaBis <- qt(1-alpha*k/2/m, df = n-2)/sqrt(n-2+qt(1-alpha*k/2/m, df = n-2)^2)
# lambdaSeries <- c(rep(lambdaBis[1], p), rep(lambdaBis, each=2))
# stock.gslope <- gslopeADMM(cov(Y), lambda = lambdaSeries)$precisionMatrix
# stock.gslope.npn <- gslopeADMM(cov(Y.npn), lambda = lambdaSeries)$precisionMatrix
# 
# plotMatrix(as.matrix(out.npn$icov[[1]]))
# plotMatrix(as.matrix(out$icov[[1]]))
# plotMatrix(stock.gslope)
# plotMatrix(stock.gslope.npn)
# plotMatrix(solve(cov(Y)))
# plotMatrix(solve(cov(Y.npn)))