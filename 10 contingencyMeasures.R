# Written by Micha³ Makowski

# Contingency measures

FP <- function(estimatedMatrix, realMatrix)
{
    # H_0: x_ij == 0
    FP <- sum(estimatedMatrix[realMatrix == 0] != 0)

    return(FP)
}

TP <- function(estimatedMatrix, realMatrix)
{
    # H_0: x_ij == 0
    TP <- sum(estimatedMatrix[realMatrix != 0] != 0)
    
    return(TP)
}

FN <- function(estimatedMatrix, realMatrix)
{
    # H_0: x_ij == 0
    FN <- sum(estimatedMatrix[realMatrix != 0] == 0)
    
    return(FN)
}

TN <- function(estimatedMatrix, realMatrix)
{
    # H_0: x_ij == 0
    TN <- sum(estimatedMatrix[realMatrix == 0] == 0)
    
    return(TN)
}

FDP <- function(estimatedMatrix, realMatrix)
{
    # H_0: x_ij == 0
    predictedPositive <- max(c(sum(estimatedMatrix != 0), 1))
    
    return(FP(estimatedMatrix, realMatrix)/predictedPositive)
}


SN <- function(estimatedMatrix, realMatrix)
{
    # H_0: x_ij == 0
    realPositive <- max(c(sum(realMatrix != 0), 1))
    
    return(TP(estimatedMatrix, realMatrix)/realPositive)
}