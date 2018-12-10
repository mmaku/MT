# Written by Micha? Makowski

# Plotting functions
require(ggplot2, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(huge, quietly = TRUE)
require(reshape2, quietly = TRUE)
require(R.utils, quietly = TRUE)
require(tikzDevice, quietly = TRUE)
require(dplyr, quietly = TRUE)

source("lambdas.R")

lambdaList <- list()
p <- 100

lambdaList$banerjee   <- lambdaSelector(p, 100, 0.05, method = "glasso")
lambdaList$holm   <- lambdaSelector(p, 100, 0.05, method = "holm")[seq(from = p+1, to = p^2, by = 2)]
lambdaList$BH   <- lambdaSelector(p, 100, 0.05, method = "BH")[seq(from = p+1, to = p^2, by = 2)]


df <- rbind(data.frame(X1 = 0:1, X2 = lambdaList$banerjee, method = "gLasso (Banerjee)"),
            data.frame(X1 = (1:length(lambdaList$BH))/length(lambdaList$BH), X2 = lambdaList$BH, 
                       method = "gSLOPE (BH)"), 
            data.frame(X1 = (1:length(lambdaList$holm))/length(lambdaList$holm), X2 = lambdaList$holm, 
                       method = "gSLOPE (Holm)"), 
            stringsAsFactors = TRUE)

ggplot(df, aes(x=X1, y=X2, color = method)) +
    geom_line() +
    geom_line(size = 1) +
    scale_color_discrete(name = "Procedure:") +
    ylim(min(lambdaList$BH),lambdaList$banerjee) +
    xlim(0,1) +
    theme_bw(base_size = 8) +
    guides(color = guide_legend(override.aes = aes(size = 2))) +
    theme(aspect.ratio = 1, 
          legend.position = "right",
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.spacing = unit(0, "cm"),
          legend.margin = margin(c(0,0,0,0)),
          plot.margin = margin(c(0,0,0,0)),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) -> myPlot


ggsave("!01 Plots/02 Other/Lambdas.png", myPlot, 
       width = 5.4, height = 5.4*myPlot$theme$aspect.ratio)

tikzTitle <- "!01 Plots/02 Other/Lambdas.tikz"

tikz(file = tikzTitle, 
     width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
plot(myPlot)
dev.off()

lines <- readLines(con = tikzTitle)
lines <- lines[-which(grepl("\\path\\[clip\\]*", x = lines, perl=F))]
lines <- lines[-which(grepl("\\path\\[use as bounding box*", x = lines, perl=F))]

linesCopy <- lines
writeLines(lines,con = tikzTitle)
