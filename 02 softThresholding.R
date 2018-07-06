# Written by Micha³ Makowski

# Soft thresholding functions definitions

softThresholding <- function(matrix,
                             threshold, 
                             penalizeDiagonal = TRUE)
{
    out <- matrix * 0
    
    out[matrix >  threshold] <- matrix[matrix >  threshold]-threshold
    out[matrix < -threshold] <- matrix[matrix < -threshold]+threshold
    
    if(!penalizeDiagonal) diag(out) <- diag(matrix)
    
    out
}

softThresholding2 <- function(matrix, 
                              threshold)
{
    pmax(matrix-threshold, 0) - pmax(-matrix-threshold, 0)
}