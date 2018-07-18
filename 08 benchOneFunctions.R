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
                     penalizeDiagonal = FALSE, 
                     truncate = TRUE, 
                     epsilon = 10e-4, 
                     mu = 1, 
                     selectCriterion = "stars")
{
    sampleCovariance <- cov(scale(data))
    
    n     <- nrow(data)
    p     <- ncol(data)
    
    banerjeeLassoLambda <- lambdaSelector(p, n, alpha, "banerjee", verbose = FALSE)
    # gLASSO
    gLassoADMM <- glassoADMM(sampleCovariance, mu = mu, lambda = banerjeeLassoLambda, 
                             penalizeDiagonal = penalizeDiagonal, truncate = truncate, absoluteEpsilon = epsilon)
    
    BHSlopeLambda <- banerjeeLassoLambda <- lambdaSelector(p, n, alpha, "BH", verbose = FALSE)
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


plotFour <- function(benchResult)
{
    benchResult <- X
    properData <- rbind(meltingMatrix((benchResult[[1]]$matrix > 0), benchResult[[1]]$name),
                        meltingMatrix((benchResult[[2]]$matrix > 0), benchResult[[2]]$name),
                        meltingMatrix((benchResult[[3]]$matrix > 0), benchResult[[3]]$name),
                        meltingMatrix((benchResult[[4]]$matrix > 0), benchResult[[4]]$name))
    
    colnames(properData) <- c("X1", "X2", colnames(properData)[3:4])
    properData$value <- c("0", "1")[1+properData$value]
    
    mainTitle <- paste0("Comparision of different precision matrix estimation methods")
    
    out <- ggplot(properData, aes(x=X1, y=X2)) + 
        geom_tile(aes(fill = factor(value))) +
        labs(title =  mainTitle, x = TeX('$X_1$'), y = TeX('$X_2$')) +
        facet_wrap(~index) +
        scale_fill_manual(values = c("lightblue", "darkblue")) +
        theme_minimal()
    out
    return(out)
}
