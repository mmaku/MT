# Written by Micha³ Makowski

require(reshape2, quietly = TRUE) # melt

source("lambdas.R")

# Melts matrix into 3-column DF
meltingMatrix <- function(matrix, 
                          index = NULL)
{
    colnames(matrix) <- seq_len(NCOL(matrix))
    rownames(matrix) <- seq_len(NROW(matrix))
    
    matrix <- melt(matrix)
    
    if(!is.null(index))
        matrix <- cbind(matrix, index)
    
    return(matrix)
}

# Simple upper triangle function
upper <- function(input)
{
    return(input[upper.tri(input, diag = FALSE)])
}

# Round data frame
round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(clusterTiming, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
}

# Adjacent matrix conversion
properAdjacent <- function(input)
{
    output <- as.matrix(input != 0)
    if(sum(diag(output)) == 0)
        diag(output) <- TRUE
    
    return(output)
}

