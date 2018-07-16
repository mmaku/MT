# Written by Micha³ Makowski

# Soft thresholding functions

softThresholding <- function(matrix,
                             threshold, 
                             penalizeDiagonal = TRUE)
{
    out <- matrix * 0
    
    out[matrix >  threshold] <- matrix[matrix >  threshold]-threshold
    out[matrix < -threshold] <- matrix[matrix < -threshold]+threshold
    
    if(!penalizeDiagonal) diag(out) <- diag(matrix)
    
    return(out)
}

softThresholding2 <- function(matrix, 
                              threshold)
{
    return(pmax(matrix-threshold, 0) - pmax(-matrix-threshold, 0))
}