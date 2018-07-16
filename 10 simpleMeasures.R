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