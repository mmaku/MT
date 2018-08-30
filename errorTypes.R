# Written by Micha³ Makowski

source("auxilaryFunctions.R") # upper
require(igraph, quietly = TRUE) # edge_connectivity, graph_from_adjacency_matrix

# Measures
# H_0: x_ij == 0
# H_1: x_ij != =

# False positive
FP <- function(estimatedMatrix, 
               adjacentMatrix)
{
    estimatedMatrix <- upper(estimatedMatrix)
    adjacentMatrix <- upper(adjacentMatrix)

    FP <- sum((estimatedMatrix != 0) & !adjacentMatrix)
    
    return(FP)
}

# True positive
TP <- function(estimatedMatrix, 
               adjacentMatrix)
{
    estimatedMatrix <- upper(estimatedMatrix)
    adjacentMatrix <- upper(adjacentMatrix)
    
    TP <- sum((estimatedMatrix != 0) & adjacentMatrix)
    
    return(TP)
}

# False negative
FN <- function(estimatedMatrix, 
               adjacentMatrix)
{
    estimatedMatrix <- upper(estimatedMatrix)
    adjacentMatrix <- upper(adjacentMatrix)
    
    FN <- sum((estimatedMatrix == 0) & adjacentMatrix)
    
    return(FN)
}

# True negative

TN <- function(estimatedMatrix, 
               adjacentMatrix)
{
    estimatedMatrix <- upper(estimatedMatrix)
    adjacentMatrix <- upper(adjacentMatrix)
    
    TN <- sum((estimatedMatrix == 0) & !adjacentMatrix)
    
    return(TN)
}

# False dicovery proportion
FDP <- function(estimatedMatrix, 
                adjacentMatrix)
{
    predictedPositive <- max(c(sum(upper(estimatedMatrix) != 0), 1))
    
    return(FP(estimatedMatrix, adjacentMatrix)/predictedPositive)
}

# Local FDP
localFDP <- function(estimatedMatrix, 
                     adjacentMatrix)
{
    adjacentMatrixIG <- graph_from_adjacency_matrix(adjacentMatrix)
    
    for (i in seq_len(NCOL(adjacentMatrix)-1)+1) 
    {
        for (j in seq_len(i-1))
        {
            adjacentMatrix[j,i] <- !isZero(edge_connectivity(adjacentMatrixIG, j, i))
        }
        
    }

    predictedPositive <- max(c(sum(upper(estimatedMatrix) != 0), 1))
    
    return(FP(estimatedMatrix, adjacentMatrix)/predictedPositive)
}


# Sensitivity
SN <- function(estimatedMatrix, 
               adjacentMatrix)
{
    realPositive <- max(c(sum(upper(adjacentMatrix)), 1))
    
    return(TP(estimatedMatrix, adjacentMatrix)/realPositive)
}

# Specificity
SP <- function(estimatedMatrix, 
               adjacentMatrix)
{
    realNegative <- max(c(sum(!upper(adjacentMatrix)), 1))
    
    return(TN(estimatedMatrix, adjacentMatrix)/realNegative)
}