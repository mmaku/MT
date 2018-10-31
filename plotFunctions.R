# Written by Micha? Makowski

# Plotting functions
require(ggplot2, quietly = TRUE)
require(latex2exp, quietly = TRUE)
require(tidyr, quietly = TRUE)


# Covariance matrix plot
plotCovarianceStructure <- function(covarianceMatrix)
{
    colnames(covarianceMatrix) <- 1:ncol(covarianceMatrix)
    rownames(covarianceMatrix) <- 1:ncol(covarianceMatrix)
    df1 <- melt(covarianceMatrix)
    df1$value <- as.factor(df1$value)
    names(df1)[1:2] <- c("X1", "X2")
    
    covPlot <- ggplot(df1, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) +
        theme(axis.ticks = element_blank(), 
              axis.text.x = element_text(angle = 330, hjust = 0, vjust = 1, colour = "grey50"),
              axis.text.y = element_text(colour = "grey50"),
              panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
              panel.background = element_rect(fill = "white", colour = "white"))
    
    return(covPlot)
}


plotMatrix <- function(matrix, 
                       title = NULL)
{
    matrix <- abs(matrix)
    colnames(matrix) <- 1:ncol(matrix)
    rownames(matrix) <- 1:ncol(matrix)
    properData <- melt(matrix)
    
    colnames(properData) <- c("X1", "X2", "value")
    
    matrixPlot <- ggplot(properData, aes(x=X1, y=X2))
    
    if(length(subset(properData,  isZero.default(value, neps=10))[[1]]) == 0)
    {
        matrixPlot <- matrixPlot +
            geom_tile(data = properData, aes(fill = value))
    } else
    {
        matrixPlot <- matrixPlot +
            geom_tile(data = subset(properData, !isZero.default(value, neps=10)), aes(fill = value)) +
            geom_tile(data = subset(properData,  isZero.default(value, neps=10)), aes(colour = "0"), 
                      linetype = 0, fill = "grey50", alpha = .5)
    }
    matrixPlot <- matrixPlot +
        labs(title = title, x = TeX('$X_1$'), y = TeX('$X_2$')) +
        scale_fill_gradient(name="Matrix\nentry\nabsolute value",
                            limits = c(.Machine$double.eps, NA)) +
        scale_colour_discrete(name=NULL) +
        guides(fill = guide_colorbar(order = 1, barwidth = 1, barheight = 10), 
               colour = guide_legend(order = 2, keywidth = 1, keyheight = 1, title.position = "bottom")) +
        theme_minimal() +
        theme(legend.spacing.y = unit(-0.3, "cm"))
    
    return(matrixPlot)
}

plotFour <- function(benchResult)
{
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

plotDiffrence <- function(estimatedMatrix, adjacentMatrix, method = "graph", graphType = NULL, n = NULL, p = NULL, alpha = NULL)
{
    properEstimated <- properAdjacent(estimatedMatrix)   
    properAdjacent <- properAdjacent(adjacentMatrix)   
    
    errorMatrix <- 1*(!properAdjacent & properEstimated) + # FP
        2*(properAdjacent & properEstimated) +             # TP
        3*(properAdjacent & !properEstimated) +            # FN
        4*(!properAdjacent & !properEstimated)             # TN
    
    properData <- meltingMatrix(errorMatrix, 
                                # NULL)
                                "Status")
    
    colnames(properData) <- c("X1", "X2", colnames(properData)[3:4])
    properData$value <- c("FP", "TP", "FN", "TN")[properData$value]
    
    x <- as.factor(properData$value)
    
    levels(x) <- c("FP", "TP", "FN", "TN")
    
    mainTitle <- paste0("Difference beetwen real & estimated matrix for ", method)
    subTitle <- paste0('graphType = ', graphType, ', n = ', n, ', p = ', p, ', alpha = ', alpha)
    
    out <- ggplot(properData, aes(x=X1, y=X2)) + 
        geom_tile(aes(fill = factor(value))) +
        # scale_fill_hue(l=50) +
        # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(2:5)]) +
        scale_fill_brewer(palette = "Spectral") +
        labs(title =  mainTitle, subtitle = subTitle, fill = "Matrix\nestimator", 
             x = TeX('$X_1$'), y = TeX('$X_2$')) +
        theme_minimal()
    
    return(out)
}

plotTwo <- function(estimatedMatrix, adjacentMatrix, method = "graph", graphType = NULL, n = NULL, p = NULL, alpha = NULL)
{
    properEstimated <- properAdjacent(estimatedMatrix)   
    properAdjacent <- properAdjacent(adjacentMatrix)   
    
    errorMatrix <- 1*(!properAdjacent & properEstimated) + # FP
        2*(properAdjacent & properEstimated) +             # TP
        3*(properAdjacent & !properEstimated) +            # FN
        4*(!properAdjacent & !properEstimated)             # TN
    
    properData <- meltingMatrix(errorMatrix, 
                                # NULL)
                                "Status")
    
    colnames(properData) <- c("X1", "X2", colnames(properData)[3:4])
    properData$value <- c("FP", "TP", "FN", "TN")[properData$value]
    
    x <- as.factor(properData$value)
    
    levels(x) <- c("FP", "TP", "FN", "TN")
    
    mainTitle <- paste0("Difference beetwen real & estimated matrix for ", method)
    subTitle <- paste0('graphType = ', graphType, ', n = ', n, ', p = ', p, ', alpha = ', alpha)
    
    out <- ggplot(properData, aes(x=X1, y=X2)) + 
        geom_tile(aes(fill = factor(value))) +
        labs(title =  mainTitle, subtitle = subTitle, x = TeX('$X_1$'), y = TeX('$X_2$')) +
        labs(fill = "Matrix\nestimator") +
        theme_minimal()
    
    return(out)
}
