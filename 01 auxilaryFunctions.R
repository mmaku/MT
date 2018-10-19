# Written by Micha³ Makowski

require(ggplot2, quietly = TRUE)
require(reshape2, quietly = TRUE)
# require(MASS, quietly = TRUE)
require(latex2exp, quietly = TRUE)

# lambdas
# meltingMatrix auxilaryFunctions
# upper auxilaryFunctions

# Covariance matrix plot
plotCovarianceStructure <- function(covarianceMatrix)
{
	colnames(covarianceMatrix) <- 1:NCOL(covarianceMatrix)
	rownames(covarianceMatrix) <- 1:NCOL(covarianceMatrix)
	df1 <- melt(covarianceMatrix)
	df1$value <- as.factor(df1$value)
	
	out <- ggplot(df1, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) +
	    theme(axis.ticks = element_blank(), 
	          axis.text.x = element_text(angle = 330, hjust = 0, vjust = 1, colour = "grey50"),
	          axis.text.y = element_text(colour = "grey50"),
	          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
	          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
	          panel.background = element_rect(fill = "white", colour = "white"))
	
	return(out)
}

# Simple matrix plot
plotMatrix <- function(matrix, 
                       title = NULL)
{
	colnames(matrix) <- 1:NCOL(matrix)
	rownames(matrix) <- 1:NCOL(matrix)
	properData <- melt(matrix)
	
	colnames(properData) <- c("X1", "X2", "value")

    out <- ggplot(properData, aes(x=X1, y=X2)) + 
        geom_tile(data = subset(properData, !(value == 0)), aes(fill = value)) +
        geom_tile(data = subset(properData,  (value == 0)), aes(colour = "0"), 
                  linetype = 0, fill = "grey50", alpha = .5) +
	    labs(title = title, x = TeX('$X_1$'), y = TeX('$X_2$')) +
        scale_fill_gradient(name="Matrix\nentry\nvalue",
                            limits = c(.Machine$double.eps, NA)) +
        scale_colour_discrete(name=NULL) +
        guides(fill = guide_colorbar(order = 1, barwidth = 1, barheight = 10), 
               colour = guide_legend(order = 2, keywidth = 1, keyheight = 1, title.position = "bottom")) +
	    theme_minimal() +
        theme(legend.spacing.y = unit(-0.3, "cm"))
    
    return(out)
}

# Bench results plot
plotFour <- function(benchResult)
{
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

    return(out)
}
