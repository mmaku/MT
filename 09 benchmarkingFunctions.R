# Written by Micha³ Makowski

# install.packages("MASS")
# install.packages("glasso")
# install.packages("huge")

require(MASS)
require(glasso)
require(huge)
require(ggplot2)
require(latex2exp)
require(reshape2)

source("01 auxilaryFunctions.R")
source("04 admmGLASSO.R")
source("07 admmGSLOPE.R")



benchOne <- function(data, 
                     alpha = .05, 
                     penalizeDiagonal = F, 
                     truncate = T, 
                     epsilon = 10e-4, 
                     mu = 1, 
                     selectCriterion = "stars")
{
    sampleCovariance <- cov(data)
    
    n     <- nrow(data)
    p     <- ncol(data)
    m     <- p*(p-1)/2
    mBanerjee <- p^2   
    banerjeeLassoLambda <- qt(1-alpha/2/mBanerjee, df = n-2)/sqrt(n-2+qt(1-alpha/2/mBanerjee, df = n-2)^2)
    
    # gLASSO
    gLassoADMM <- glassoADMM(sampleCovariance, mu = mu, lambda = banerjeeLassoLambda, 
                             penalizeDiagonal = penalizeDiagonal, truncate = truncate, absoluteEpsilon = epsilon)
    
    k = 1:m
    BHlambda <- qt(1-alpha*k/2/m, df = n-2)/sqrt(n-2+qt(1-alpha*k/2/m, df = n-2)^2)
    BHSlopeLambda <- c(rep(BHlambda[1], p), rep(BHlambda, each=2))
    # gSLOPE
    gSlopeADMM <- gslopeADMM(sampleCovariance, mu = mu, lambda = BHSlopeLambda, 
                             penalizeDiagonal = penalizeDiagonal, truncate = truncate, absoluteEpsilon = epsilon)
    
    # gLASSOstars
    gLassoHuge <- huge(data, nlambda=30, method = "glasso")
    lassoSelect <- huge.select(gLassoHuge, criterion = selectCriterion) # Select the graph using StARS
    gLassoStARS <- glassoADMM(sampleCovariance, mu = mu, lambda = lassoSelect$opt.lambda, 
                              penalizeDiagonal = penalizeDiagonal, truncate = truncate, absoluteEpsilon = epsilon)
    
    return(list(gLassoADMM = list(lassoParameter = gLassoADMM$lambda, 
                                  precisionMatrix = gLassoADMM$precisionMatrix, 
                                  residuals = gLassoADMM$residuals, 
                                  iterations = gLassoADMM$iterations),
                gSlopeADMM = list(lassoParameter = gSlopeADMM$lambda, 
                                  precisionMatrix = gSlopeADMM$precisionMatrix, 
                                  residuals = gSlopeADMM$residuals, 
                                  iterations = gSlopeADMM$iterations),
                gLassoStARS = list(lassoParameter = gLassoStARS$lambda, 
                                   precisionMatrix = gLassoStARS$precisionMatrix, 
                                   residuals = gLassoStARS$residuals, 
                                   iterations = gLassoStARS$iterations),
                sampleCovariance = sampleCovariance,
                lagrangianParameter = mu,
                diagonalPenalization = penalizeDiagonal,
                truncated = truncate,
                epsilon = epsilon,
                alpha = alpha,
                selectMethod = selectCriterion))
}


strictGraphicalFDR <- function(estimatedMatrix, realMatrix)
{
    # H_0: x_ij == 0
    falsePositive <- sum(estimatedMatrix[realMatrix == 0] != 0)
    predictedPositive <- sum(estimatedMatrix != 0)
    
    return(falsePositive/predictedPositive)
}

x <- matrix(c(1,2,3,0,0,1,2,3,4,0,0,4,2,2,4,5), nrow = 4, ncol = 4)
sum(t(x)[x == 0] != 0)


meltingMatrix <- function(matrix, index = NULL)
{
    colnames(matrix) <- 1:ncol(matrix)
    rownames(matrix) <- 1:ncol(matrix)
    
    matrix <- reshape2::melt(matrix)
    
    if(!is.null(index))
        matrix <- cbind(matrix, index)
    
    return(matrix)
}
    
plotBenchOne <- function(benchResult)
{
    properData <- rbind(meltingMatrix(benchResult[[1]]$precisionMatrix, "gLASSO"),
                        meltingMatrix(benchResult[[2]]$precisionMatrix, "gSLOPE"),
                        meltingMatrix(benchResult[[3]]$precisionMatrix, "gLASSO_StARS"),
                        meltingMatrix(benchResult$sampleCovariance, "sampleCovariance"))

    colnames(properData) <- c("X1", "X2", colnames(properData)[3:4])
    
    mainTitle <- paste0("Comparision of different precision matrix estimation methods")
    
    subTitle <- paste0('alpha = ', benchResult$alpha, 
                        ', diagonal penalization: ', c('No', 'Yes')[1+benchResult$diagonalPenalization],
                        ', values below ', benchResult$epsilon, ' truncated: ', c('No.', 'Yes.')[1+benchResult$truncated])
    
    properData
    
    out <- ggplot(properData, aes(x=X1, y=X2))
    
    if(benchResult$truncated)
    {
        out = out + geom_tile(data = subset(properData, !is.zero(value)), aes(fill = value)) + 
            geom_tile(data = subset(properData,  is.zero(value)), aes(colour = "0"), linetype = 0, fill = "grey50", alpha = .5)
    }
    else
    {
        out = out + geom_tile(aes(fill = value))
    }
    
    out = out + labs(title =  mainTitle, subtitle = subTitle, x = TeX('$X_1$'), y = TeX('$X_2$')) +
        facet_wrap(~index) +
        scale_fill_gradient(name="Matrix\nentry\nvalue",
                            limits = c(.Machine$double.eps, NA)) +
        scale_colour_discrete(name=NULL) +
        guides(fill = guide_colorbar(order = 1, barwidth = 1, barheight = 10), 
               colour = guide_legend(order = 2, keywidth = 1, keyheight = 1, title.position = "bottom")) +
        theme_minimal() +
        theme(legend.spacing.y = unit(-0.3, "cm"))

    return(out)
}
