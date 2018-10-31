# Written by Micha³ Makowski

require(dplyr, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(tikzDevice, quietly = TRUE)
require(xtable, quietly = TRUE)


load("./!02 Data/01 Binded/01 AllOne/ROC.RData")

results %>%
    select(-c(procedure, multi)) %>%
    as_tibble() -> finalResults

colnames(finalResults)[5] <- "graph.u"
colnames(finalResults)[3] <- "graphType"

gt <- unique(finalResults$graphType)[1]
u <- unique(finalResults$graph.u)[1]

for(gt in unique(finalResults$graphType))
{
    for(u in unique(finalResults$graph.u))
    {
        finalResults %>%
            filter(graphType == gt & graph.u == u) %>%
            mutate(met = recode_factor(met,
                                             `banerjee` = "gLasso (Banerjee)", 
                                             `BH` = "gSLOPE (BH)",
                                             `holm` = "gSLOPE (Holm)")) %>%
            ggplot(aes(x = 1-SP, y = SN, color = met)) +
            geom_line(size = 1)+ 
            ylim(c(0,1)) +
            xlim(c(0,1)) +
            labs(title = "ROC curve",
                 subtitle = paste("Setup: p = $100", 
                                  ifelse(gt != "scale-free", 
                                         paste("$, number of components = $10$,"), 
                                         "$,"), 
                                  "scaled parameter $\\alpha$ = $0.05$, graph type =", gt,
                                  ifelse(gt == "cluster", 
                                         paste(", $\\prob(x_{ij}\\neq0)=1$,"), 
                                         ","),
                                  "SNR = $",round(1/(1-u),2),"$."),
                 y = "Sensitivity",
                 x = "1-Specificity") +
            labs(subtitle = paste0("Setup: scaled $\\alpha$ = $", a, "$, $",  p,
                                   ifelse(gt != "scale-free", 
                                          paste0("$ variables, $", g, "$ components, "), 
                                          "$ variables, "),
                                   gt, " graph",
                                   ifelse(gt == "cluster", 
                                          paste0(", $prob(x_{ij}\\neq0)=", prob, "$."), 
                                          ".")),
                 y = "Value",
                 x = "Sample size n") +            
            scale_color_discrete(name = "Procedure:") +
            # theme_bw() +
            theme_bw(base_size = 8) +
            guides(color = guide_legend(override.aes = aes(size = 2)),
                   shape = guide_legend(override.aes = aes(size = 2))) +
            theme(aspect.ratio = 8/8, 
                  plot.margin = margin(c(0,0,0,0)),
                  legend.margin = margin(c(0,10,0,0)))  -> myPlot
        
        ggsave(paste0("!01 Plots/01 Results/04 ROC/ROCprob1_", gt, "_", u, ".png"), myPlot, 
               width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
        dev.off()
        
        tikz(file = paste0("!01 Plots/01 Results/04 ROC/ROCprob1_", gt, "_", u, ".tikz"), 
             width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
        plot(myPlot)
        dev.off()
        
        # print(xtable(myResults,
        #              caption = "Power and FDR of each procedure",
        #              auto = TRUE),
        #       booktabs = TRUE)
    }
}
