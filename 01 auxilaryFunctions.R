# Written by Micha³ Makowski


# TODO
# Penalize diagonal in lambdas functions

require(ggplot2)
require(reshape)
require(MASS)
require(latex2exp)

plotCovarianceStructure <- function(covarianceMatrix)
{
	colnames(covarianceMatrix) <- 1:ncol(covarianceMatrix)
	rownames(covarianceMatrix) <- 1:ncol(covarianceMatrix)
	df1 <- melt(covarianceMatrix)
	df1$value <- as.factor(df1$value)
	
	covPlot <- ggplot(df1, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) +
	    theme(axis.ticks = element_blank(), 
	          axis.text.x = element_text(angle = 330, hjust = 0, vjust = 1, colour = "grey50"),
	          axis.text.y = element_text(colour = "grey50"),
	          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
	          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
	          panel.background = element_rect(fill = "white", colour = "white"))
	
	return(covPlot)
}

is.zero <- function(x)
{
    return(x == 0)
}

plotMatrix <- function(matrix, 
                       title = NULL)
{
	colnames(matrix) <- 1:ncol(matrix)
	rownames(matrix) <- 1:ncol(matrix)
	properData <- melt(matrix)
	
	colnames(properData) <- c("X1", "X2", "value")

    matrixPlot <- ggplot(properData, aes(x=X1, y=X2)) + 
        geom_tile(data = subset(properData, !is.zero(value)), aes(fill = value)) +
        geom_tile(data = subset(properData,  is.zero(value)), aes(colour = "0"), 
                  linetype = 0, fill = "grey50", alpha = .5) +
	    labs(title = title, x = TeX('$X_1$'), y = TeX('$X_2$')) +
        scale_fill_gradient(name="Matrix\nentry\nvalue",
                            limits = c(.Machine$double.eps, NA)) +
        scale_colour_discrete(name=NULL) +
        guides(fill = guide_colorbar(order = 1, barwidth = 1, barheight = 10), 
               colour = guide_legend(order = 2, keywidth = 1, keyheight = 1, title.position = "bottom")) +
	    theme_minimal() +
        theme(legend.spacing.y = unit(-0.3, "cm"))
    
    return(matrixPlot)
}

# Banerjee for gLASSO

lambdaBanerjee <- function(input, n, alpha = 0.05, verbose = TRUE)
{
    p <- ncol(input)
    
    if(!is.matrix(input))
    {
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

    pBanerjee <- p^2
    fraction <- qt(1-alpha/2/pBanerjee, df = n-2)/sqrt(n-2+qt(1-alpha/2/pBanerjee, df = n-2)^2)

    return(twoLargestProd*fraction)
}

# Bonferonni for gLASSO (not sure if done properly - problem with sigma i.e. twoLargestProd if data not scaled)

lambdaBonf <- function(input, n, alpha = 0.05, verbose = TRUE)
{
    p <- ncol(input)
    
    if(!is.matrix(input))
    {
        p <- input[1]
        twoLargestProd <- 1
        
    } else if(!isSymmetric(input))
    {
        if(verbose) cat("The input is identified as the data matrix.\n")
        
        input <- cov(scale(input))
        twoLargestProd <- 1
    } else
    {
        if(verbose) cat("The input is identified as the covariance matrix.\n")
        
        twoLargestProd <- prod(-sort(-diag(input), partial = 2)[1:2]) # In case data wasn't scaled
    }
    
    qBonf <- qnorm(1-alpha/2/p)
    
    return(twoLargestProd*qBonf)
}

# Banerjee SLOPE modification

lambdaBS <- function(input, n, alpha = 0.05, verbose = TRUE)
{
    p <- ncol(input)
    
    if(!is.matrix(input))
    {
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
    
    pBS <- p*(p-1)/2
    k <- 1:pBS
    fractionSeq <- qt(1-alpha*k/2/pBS, df = n-2)/sqrt(n-2+qt(1-alpha*k/2/pBS, df = n-2)^2)
    fractionSeq <- c(rep(fractionSeq[1], p), rep(fractionSeq, each=2))
    
    return(twoLargestProd*fractionSeq)
}

# BH for SLOPE (not sure if done properly - problem with sigma i.e. twoLargestProd if data not scaled)

lambdaBH <- function(input, n, alpha = 0.05, verbose = TRUE)
{
    p <- ncol(input)
    
    if(!is.matrix(input))
    {
        p <- input[1]
        twoLargestProd <- 1
        
    } else if(!isSymmetric(input))
    {
        if(verbose) cat("The input is identified as the data matrix.\n")
        
        input <- cov(scale(input))
        twoLargestProd <- 1
    } else
    {
        if(verbose) cat("The input is identified as the covariance matrix.\n")
        
        twoLargestProd <- prod(-sort(-diag(input), partial = 2)[1:2]) # In case data wasn't scaled
    }
    
    pBH <- p*(p-1)/2
    k <- 1:pBH
    fractionSeq <- qnorm(1-alpha*k/2/pBH)
    fractionSeq <- c(rep(fractionSeq[1], p), rep(fractionSeq, each=2))
    
    return(twoLargestProd*fractionSeq)
}

lambdaSelector <- function(input, n, alpha = 0.05, method = "banerjee", verbose = TRUE)
{
    out = switch(method,
                 banerjee = lambdaBanerjee(input, n, alpha, verbose),
                 bonf = lambdaBonf(input, n, alpha, verbose),
                 BS = lambdaBS(input, n, alpha, verbose),
                 BH = lambdaBH(input, n, alpha, verbose))
                 
    return(out)
}

localFDR <- function(matrix, 
                     res)
{
	p=ncol(matrix)

	
		sapply(1:p, function(i) sum(matrix[,i]==0 & res[,i]!=0)/max(1, sum(res[,i]!=0)-1))
}

dataSimulator <- function(n = 100, 
                          SNR = 1, 
                          K = 10, 
                          numb.vars = 30, 
                          max.dim = 2, 
                          equal.dims = TRUE)
{         
    sigma <- 1/SNR
    if (equal.dims)
        dims <- rep(max.dim, K)
    else dims <- sample(1:max.dim, K, replace = T)
    X <- NULL
    Y <- NULL
    s <- NULL
    factors <- NULL
    for (j in 1:K) 
    {
        Z <- qr.Q(qr(replicate(dims[j], rnorm(n, 0, 1))))
        coeff <- matrix(runif(dims[j] * numb.vars, 0.1, 1) *
                            sign(runif(dims[j] * numb.vars, -1, 1)), nrow = dims[j])
        SIGNAL <- Z %*% coeff
        SIGNAL <- scale(SIGNAL)
        Y <- cbind(Y, SIGNAL)
        factors <- cbind(factors, Z)
        X <- cbind(X, SIGNAL + replicate(numb.vars, rnorm(n, 0, sigma)))
        s <- c(s, rep(j, numb.vars))
    }
    
    return(list(X = X, 
                signals = Y, 
                factors = factors, 
                dims = dims, 
                s = s))
}