# Written by Micha? Makowski

# Plotting functions
require(ggplot2, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(huge, quietly = TRUE)
require(reshape2, quietly = TRUE)
require(R.utils, quietly = TRUE)
require(tikzDevice, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(colorRamps, quietly = TRUE)

modelList <- list()

# modelList$one   <- list(v = 0.2, u = 0.8, prob = 0.5)
modelList$two   <- list(v = 0.3, u = 0.7, prob = 0.5)
# modelList$three <- list(v = 0.4, u = 0.6, prob = 0.5)
# modelList$four  <- list(v = 0.6, u = 0.4, prob = 0.5)
modelList$five  <- list(v = 0.7, u = 0.3, prob = 0.5)
# modelList$six   <- list(v = 0.8, u = 0.2, prob = 0.5)
modelList$seven <- list(v = -1, u = 1, prob = 0.5)
# modelList$eight <- list(v = -0.5, u = 1, prob = 0.5)

df <- data.frame(X1 = numeric(), X2 = numeric(), value = numeric(), 
                 MR = character(), type = character(), stringsAsFactors = TRUE)

m <- modelList[[1]]
for(m in modelList) 
{
    MR <- round((m$v + m$u)/m$v, digits = 2)
    
    model <- doCall("huge.generator", 
                    n = 30, d = 30, g=3, graph = "cluster", verbose = FALSE, vis = FALSE,
                    args = m)
    
    df <- rbind(df,
                cbind(melt(model$sigmahat), MR, type = "Sample covariance matrix"),
                cbind(melt(model$sigma), MR, type = "True covariance matrix"),
                cbind(melt(model$omega), MR, type = "Precision matrix"))
}

mutate(df, MRtext = paste0("$SNR=", MR, "$")) %>%
    dplyr::select(-MR) -> df

# matrixPlot <- 
ggplot(df, aes(x=Var1, y=Var2, fill = value)) +
    geom_tile() +
    # geom_tile(data = subset(df, !isZero.default(value, neps=10)), aes(fill = exp(value))) +
    # geom_tile(data = subset(df,  isZero.default(value, neps=10)), aes(color = '0'),
              # linetype = 0, fill = "grey70", alpha = 1) +
    # labs(title = title, x = "$X_1$", y = "$X_2$") +
    # scale_fill_gradient(name="Matrix\nentry\n value",
                        # limits = c(.Machine$double.eps, NA)) +
    # scale_colour_discrete(name=NULL) +
    # guides(fill = guide_colorbar(order = 1, barwidth = 1, barheight = 10), 
    #        colour = guide_legend(order = 2, keywidth = 1, keyheight = 1, title.position = "bottom")) +
    facet_grid(type ~ MRtext) +
    labs(subtitle = "",
         y = "",
         x = "") +
    scale_fill_gradientn(colors = matlab.like(275)[(1:6)*50-25]) +
    theme_bw(base_size = 8) +
    # theme(aspect.ratio = 8/16, 
    #       plot.margin = margin(c(0,0,0,0)),
    #       legend.margin = margin(c(0,10,0,0))) -> myPlot
    theme(aspect.ratio = 12/16, 
          plot.margin = margin(c(0,0,0,0)),
          legend.position="none") -> myPlot


ggsave("!01 Plots/02 Other/MR_cluster.png", myPlot, 
       width = 5.4, height = 5.4*myPlot$theme$aspect.ratio)

tikzTitle <- "!01 Plots/02 Other/MR_cluster.tikz"

tikz(file = tikzTitle, 
     width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
plot(myPlot)
dev.off()

lines <- readLines(con = tikzTitle)
lines <- lines[-which(grepl("\\path\\[clip\\]*", x = lines, perl=F))]
lines <- lines[-which(grepl("\\path\\[use as bounding box*", x = lines, perl=F))]
lines <- gsub(pattern = "SNR", replace = "\\SNR", x = lines, fixed = TRUE)

for(lineNo in seq_along(lines))
{
    while(substr(lines[lineNo], 39, 47) == "rectangle" &
          substr(lines[lineNo], 50, 55) == substr(lines[lineNo+2], 24, 29) &
          substr(lines[lineNo], 31, 36) == substr(lines[lineNo+2], 31, 36) &
          substr(lines[lineNo], 57, 62) == substr(lines[lineNo+2], 57, 62))
    {
        substr(lines[lineNo], 50, 55) <- substr(lines[lineNo+2], 50, 55)
        lines <- lines[-(lineNo + 1:2)]
    }
    if(lineNo > length(lines)-3)
        break
}

writeLines(lines,con = tikzTitle)
#
# for(lineNo in seq_along(lines))
# {
#     if(substr(lines[lineNo], 39, 47) == "rectangle")
#     {
#         for(lineNo2 in (lineNo+1):length(lines))
#         {
#             if(substr(lines[lineNo2], 39, 47) == "rectangle" &
#                substr(lines[lineNo2-2], 30, 41) == substr(lines[lineNo-2], 30, 41) &
#                substr(lines[lineNo], 24, 29) == substr(lines[lineNo2], 24, 29) &
#                substr(lines[lineNo], 50, 55) == substr(lines[lineNo2], 50, 55) &
#                substr(lines[lineNo], 57, 62) == substr(lines[lineNo2], 31, 36))
#             {
#                 substr(lines[lineNo], 57, 62) <- substr(lines[lineNo2], 57, 62)
#                 lines <- lines[-(lineNo2 - 0:2)]
#                 lineNo2 = lineNo2-4
#             }
#             if(lineNo2 >= length(lines))
#                 break
#         }
#     }
#
#     if(lineNo > length(lines)-3)
#         break
# }
#
# tikzTitle <- "!01 Plots/02 Other/MR3.tikz"
# writeLines(lines,con = tikzTitle)
