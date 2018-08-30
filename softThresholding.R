# Written by Micha³ Makowski

# Soft thresholding function
softThresholding <- function(matrix,
                             threshold, 
                             penalizeDiagonal = TRUE)
{
    output <- matrix * 0
    
    output[matrix >  threshold] <- matrix[matrix >  threshold]-threshold
    output[matrix < -threshold] <- matrix[matrix < -threshold]+threshold
    
    if(!penalizeDiagonal) diag(output) <- diag(matrix)
    
    return(output)
}