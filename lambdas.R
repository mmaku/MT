# Written by Micha³ Makowski

# Functions for lambda regularizer 
# Banerjee for gLasso
# BH and Holm for gSLOPE

# Lambda selector/wrapper for diffrent types of procedures
lambdaSelector <- function(input, n, alpha = 0.05, method = "banerjee", verbose = TRUE)
{
    p <- ncol(input)
    
    if(alpha != 0)
    {
        if(!is.matrix(input))
        {
            if(verbose) cat("The input is identified as a dimension.\n")
            if(!is.numeric(input[1])) stop(paste("The input must be numeric, but is", typeof(input[1])))

            p <- input[1]
            twoLargestProd <- 1
            
        } else if(!isSymmetric(input))
        {
            if(verbose) cat("The input is identified as the data matrix.\n")
            
            n <- nrow(input)
            input <- cov(scale(input))
            twoLargestProd <- 1
        } else
        { 
            if(verbose) cat("The input is identified as the covariance matrix.\n")
            
            twoLargestProd <- prod(-sort(-diag(input), partial = 2)[1:2]) # In case data wasn't scaled
        }
        
        out <- switch(method,
                      glasso = lambdaGLASSO(p, n, alpha, twoLargestProd),
                      banerjee = lambdaBanerjee(p, n, alpha, twoLargestProd),
                      BH = lambdaBH(p, n, alpha, twoLargestProd),
                      holm = lambdaHolm(p, n, alpha, twoLargestProd))
    } else
    {
        out <- 0
    }
    
    return(out)
}

# gLASSO (not sure if done properly - problem with sigma i.e. twoLargestProd if data not scaled)
lambdaGLASSO <- function(p, n, alpha = 0.05, twoLargestProd = 1)
{
    pGLASSO <- p*(p-1)/2
    fraction <- qt(1-alpha/2/pGLASSO, df = n-2)/sqrt(n-2+qt(1-alpha/2/pGLASSO, df = n-2)^2)
    
    return(twoLargestProd*fraction)
}

# Banerjee for gLASSO (not sure if done properly - problem with sigma i.e. twoLargestProd if data not scaled)
lambdaBanerjee <- function(p, n, alpha = 0.05, twoLargestProd = 1)
{
    pBanerjee <- p^2
    fraction <- qt(1-alpha/2/pBanerjee, df = n-2)/sqrt(n-2+qt(1-alpha/2/pBanerjee, df = n-2)^2)
    
    return(twoLargestProd*fraction)
}

# BH for SLOPE (not sure if done properly - problem with sigma i.e. twoLargestProd if data not scaled)
lambdaBH <- function(p, n, alpha = 0.05, twoLargestProd = 1)
{
    pBH <- p*(p-1)/2
    k <- 1:pBH
    fractionSeq <- qt(1-alpha*k/2/pBH, df = n-2)/sqrt(n-2+qt(1-alpha*k/2/pBH, df = n-2)^2)
    fractionSeq <- c(rep(fractionSeq[1], p), rep(fractionSeq, each=2))
    
    return(twoLargestProd*fractionSeq)
}

# Holm for SLOPE (not sure if done properly - problem with sigma i.e. twoLargestProd if data not scaled)
lambdaHolm <- function(p, n, alpha = 0.05, twoLargestProd = 1)
{
    pHolm <- p*(p-1)/2
    k <- 1:pHolm
    fractionSeq <- qt(1-alpha/2/(pHolm+1-k), df = n-2)/sqrt(n-2+qt(1-alpha/2/(pHolm+1-k), df = n-2)^2)
    fractionSeq <- c(rep(fractionSeq[1], p), rep(fractionSeq, each=2))
    
    return(twoLargestProd*fractionSeq)
}
