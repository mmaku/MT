# Written by Micha³ Makowski


# TODO
# Lambdas

require(ggplot2, quietly = TRUE)
require(reshape, quietly = TRUE)
require(MASS, quietly = TRUE)
require(latex2exp, quietly = TRUE)


meltingMatrix <- function(matrix, 
                          index = NULL)
{
    colnames(matrix) <- 1:ncol(matrix)
    rownames(matrix) <- 1:ncol(matrix)
    
    matrix <- reshape2::melt(matrix)
    
    if(!is.null(index))
        matrix <- cbind(matrix, index)
    
    return(matrix)
}



# Simple upper triangle function

upper <- function(input, diagonal = FALSE)
{
    return(input[upper.tri(input, diag = diagonal)])
}

# Covariance matrix plot

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

# Simple zero check 

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

# Banerjee for gLASSO (not sure if done properly - problem with sigma i.e. twoLargestProd if data not scaled)

lambdaBanerjee <- function(p, n, alpha = 0.05, twoLargestProd = 1)
{
    pBanerjee <- p^2
    fraction <- qt(1-alpha/2/pBanerjee, df = n-2)/sqrt(n-2+qt(1-alpha/2/pBanerjee, df = n-2)^2)

    return(twoLargestProd*fraction)
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

# BH for SLOPE (not sure if done properly - problem with sigma i.e. twoLargestProd if data not scaled)

lambdaBH <- function(p, n, alpha = 0.05, twoLargestProd = 1)
{
    pBH <- p*(p-1)/2
    k <- 1:pBH
    fractionSeq <- qt(1-alpha*k/2/pBH, df = n-2)/sqrt(n-2+qt(1-alpha*k/2/pBH, df = n-2)^2)
    fractionSeq <- c(rep(fractionSeq[1], p), rep(fractionSeq, each=2))
    
    return(twoLargestProd*fractionSeq)
}

# Lambda selector for diffrent types of procedures

lambdaSelector <- function(input, n, alpha = 0.05, method = "banerjee", verbose = TRUE)
{
    p <- ncol(input)
    
    if(!is.matrix(input))
    {
        if(verbose) cat("The input is identified as a dimension.\n")
        
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
    
    out = switch(method,
                 banerjee = lambdaBanerjee(p, n, alpha, twoLargestProd),
                 BH = lambdaBH(p, n, alpha, twoLargestProd),
                 holm = lambdaHolm(p, n, alpha, twoLargestProd))
    
    return(out)
}

# TODO

localFDR <- function(matrix, 
                     res)
{
	p=ncol(matrix)

	
		sapply(1:p, function(i) sum(matrix[,i]==0 & res[,i]!=0)/max(1, sum(res[,i]!=0)-1))
}

# Data simulator based on SNR concept

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


plotFour <- function(benchResult)
{
    benchResult <- X
    properData <- rbind(meltingMatrix((benchResult[[1]]$matrix > 0), benchResult[[1]]$name),
                        meltingMatrix((benchResult[[2]]$matrix > 0), benchResult[[2]]$name),
                        meltingMatrix((benchResult[[3]]$matrix > 0), benchResult[[3]]$name),
                        meltingMatrix((benchResult[[4]]$matrix > 0), benchResult[[4]]$name))
    
    colnames(properData) <- c("X1", "X2", colnames(properData)[3:4])
    properData$value <- c("0", "1")[1+properData$value]
    
    mainTitle <- paste0("Comparision of different precision matrix estimation methods")
    
    out <- ggplot(properData, aes(x=X1, y=X2)) + 
        geom_tile(aes(fill = factor(value))) +
        labs(title =  mainTitle, x = TeX('$X_1$'), y = TeX('$X_2$')) +
        facet_wrap(~index) +
        scale_fill_manual(values = c("lightblue", "darkblue")) +
        theme_minimal()
    out
    return(out)
}
