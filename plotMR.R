# Written by Micha? Makowski

# Plotting functions
require(ggplot2, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(huge, quietly = TRUE)
require(reshape2, quietly = TRUE)
require(R.utils, quietly = TRUE)
require(dplyr, quietly = TRUE)

modelList <- list()

modelList$three <- list(v = 0.3, u = 0.7, prob = 0.5, g = 10)
modelList$one <- list(v = 0.7, u = 0.3, prob = 0.5, g = 10)
modelList$zero <- list(v = -1, u = 1, prob = 0.5, g = 10)

df <- data.frame(X1 = numeric(), X2 = numeric(), value = numeric(), 
                 MR = character(), type = character(), stringsAsFactors = TRUE)

m <- modelList[[1]]
for(m in modelList) 
{
    MR <- round((m$v + m$u)/m$v, digits = 2)
    
    model <- doCall("huge.generator", 
                    n = 100, d = 100, graph = "cluster", verbose = FALSE, vis = FALSE,
                    args = m)
    
    df <- rbind(df,
                cbind(melt(model$sigma), MR, type = "cov"),
                cbind(melt(model$omega), MR, type = "prec"))
}

# matrixPlot <- 
ggplot(df, aes(x=Var1, y=Var2)) +
    geom_tile(data = subset(properData, !isZero.default(value, neps=10)), aes(fill = value)) +
    geom_tile(data = subset(properData,  isZero.default(value, neps=10)), aes(colour = "0"), 
              linetype = 0, fill = "grey50", alpha = .5) +
    labs(title = title, x = "$X_1$", y = "$X_2$") +
    scale_fill_gradient(name="Matrix\nentry\n value",
                        limits = c(.Machine$double.eps, NA)) +
    scale_colour_discrete(name=NULL) +
    guides(fill = guide_colorbar(order = 1, barwidth = 1, barheight = 10), 
           colour = guide_legend(order = 2, keywidth = 1, keyheight = 1, title.position = "bottom")) +
    theme_minimal() +
    theme(legend.spacing.y = unit(-0.3, "cm")) -> plot
    
plot
