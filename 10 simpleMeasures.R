# Written by Micha³ Makowski

# Measures
# H_0: x_ij == 0
# H_1: x_ij != =

source("01 auxilaryFunctions.R")

properAdjacent <- function(input)
{
    output <- as.matrix(input != 0)
    if(sum(diag(output)) == 0)
        diag(output) <- TRUE
    
    return(output)
}

FP <- function(estimatedMatrix, 
               adjacentMatrix,
               partial = FALSE)
{
    if(partial)
    {
        estimatedMatrix <- upper(estimatedMatrix)
        adjacentMatrix <- upper(adjacentMatrix)
    }
    
    FP <- sum(estimatedMatrix[!adjacentMatrix] != 0)

    return(FP)
}

TP <- function(estimatedMatrix, 
               adjacentMatrix,
               partial = FALSE)
{
    if(partial)
    {
        estimatedMatrix <- upper(estimatedMatrix)
        adjacentMatrix <- upper(adjacentMatrix)
    }
    
    TP <- sum(estimatedMatrix[adjacentMatrix] != 0)
    
    return(TP)
}

FN <- function(estimatedMatrix, 
               adjacentMatrix,
               partial = FALSE)
{
    if(partial)
    {
        estimatedMatrix <- upper(estimatedMatrix)
        adjacentMatrix <- upper(adjacentMatrix)
    }
    
    FN <- sum(estimatedMatrix[adjacentMatrix] == 0)
    
    return(FN)
}

TN <- function(estimatedMatrix, 
               adjacentMatrix,
               partial = FALSE)
{
    if(partial)
    {
        estimatedMatrix <- upper(estimatedMatrix)
        adjacentMatrix <- upper(adjacentMatrix)
    }
    
    TN <- sum(estimatedMatrix[!adjacentMatrix] == 0)
    
    return(TN)
}

FDP <- function(estimatedMatrix, 
                adjacentMatrix,
                partial = FALSE)
{
    if(partial)
    {
        estimatedMatrix <- upper(estimatedMatrix)
        adjacentMatrix <- upper(adjacentMatrix)
    }
    
    predictedPositive <- max(c(sum(estimatedMatrix != 0), 1))
    
    return(FP(estimatedMatrix, adjacentMatrix)/predictedPositive)
}

SN <- function(estimatedMatrix, 
               adjacentMatrix,
               partial = FALSE)
{
    if(partial)
    {
        estimatedMatrix <- upper(estimatedMatrix)
        adjacentMatrix <- upper(adjacentMatrix)
    }
    
    realPositive <- max(c(sum(adjacentMatrix), 1))
    
    return(TP(estimatedMatrix, adjacentMatrix)/realPositive)
}

SP <- function(estimatedMatrix, 
               adjacentMatrix,
               partial = FALSE)
{
    if(partial)
    {
        estimatedMatrix <- upper(estimatedMatrix)
        adjacentMatrix <- upper(adjacentMatrix)
    }
    
    realNegative <- max(c(sum(!adjacentMatrix), 1))
    
    return(TN(estimatedMatrix, adjacentMatrix)/realNegative)
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
        labs(title =  mainTitle, subtitle = subTitle, x = TeX('$X_1$'), y = TeX('$X_2$')) +
        labs(fill = "Matrix\nestimator") +
        # scale_y_discrete() +
        theme_minimal()

    return(out)
}
