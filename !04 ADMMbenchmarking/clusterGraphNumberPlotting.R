# Written by Micha³ Makowski

require(dplyr, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(tikzDevice, quietly = TRUE)
require(xtable, quietly = TRUE)


load("./!02 Data/01 Binded/01 AllOne/AllOneSimulations_Num@18_10_21_13_06#384.RData")

finalResults %>%
    select(-c(localFDR, penalizeDiagonal, iterations, SP)) %>%
    filter(procedure != "gLASSO") %>%
    as_tibble() -> finalResults1

load("./!02 Data/01 Binded/01 AllOne/AllOneSimulations_Prob@18_10_21_17_24#384.RData")

finalResults %>%
    select(-c(localFDR, penalizeDiagonal, iterations, SP)) %>%
    filter(procedure != "gLASSO" & graph.prob == 0.5) %>%
    as_tibble() -> finalResults2

load("./!02 Data/01 Binded/01 AllOne/AllOneSimulations_Corr@18_10_21_11_23#576.RData")

finalResults %>%
    select(-c(localFDR, penalizeDiagonal, iterations, SP)) %>%
    filter(procedure != "gLASSO" & alpha != 0.15 & graphType == "hub" & ((graph.v == 0.7)|(graph.v == 0.3))) %>%
    as_tibble() -> finalResults3


finalResults <- rbind(finalResults1, finalResults2, finalResults3)

a <- unique(finalResults$alpha)[1]
gt <- unique(finalResults$graphType)[1]

for(gt in unique(finalResults$graphType))
{
    for(a in unique(finalResults$alpha))
    {
        p = unique(finalResults$p)
        prob = unique(finalResults$graph.prob)
        
        colLabels <- sapply(unique(finalResults$graph.g) , function(x) paste0("\\#[conectivity components] = ", x))
        names(colLabels) <- unique(finalResults$graph.g)
        
        rowLabels <- sapply(unique(finalResults$graph.v) , function(x) paste0("$SNR = ", round(1/x, 2),"$"))
        names(rowLabels) <- unique(finalResults$graph.v)
        
        finalResults %>%
            filter(graphType == gt & alpha == a) %>%
            mutate(procedure = recode_factor(procedure,
                                             `banerjee.gLASSO` = "gLasso (Banerjee)", 
                                             `BH.gSLOPE` = "gSLOPE (BH)",
                                                 `holm.gSLOPE` = "gSLOPE (Holm)")) -> myResults
        
        gather(myResults, FDR:Power, key = "metric", value = "value") %>%
            ggplot(aes(x = factor(n), y = value, color = procedure, shape = metric)) +
            geom_jitter(width = 0.2, height = 0, size = 1.5) +
            geom_hline(aes(yintercept = a)) +
            ylim(c(0,1)) +
            facet_grid(cols = vars(graph.g), rows = vars(graph.v), 
                       labeller = labeller(.rows =  rowLabels,
                                           .cols = colLabels), 
                       scales = "fixed") +
            labs(title = "Power and FDR comparison",
                 subtitle = paste("Setup: p = $", p, 
                                  "$, parameter $\\alpha$ = $", a, 
                                  "$, graph type =", gt,
                                  ifelse(gt == "cluster", 
                                         paste(", $\\P(x_{ij}\\neq0)=", prob, "$."), 
                                         ".")),
                 y = "Value",
                 x = "Sample size n") +
            scale_color_discrete(name = "Procedure:") +
            scale_shape_discrete(name = "Measure:") +
            # theme_bw() +
            theme_bw(base_size = 7) +
            guides(color = guide_legend(override.aes = aes(size = 2)),
                   shape = guide_legend(override.aes = aes(size = 2))) +
            theme(aspect.ratio = 10/16, 
                  plot.margin = margin(c(0,0,0,0)),
                  legend.margin = margin(c(0,10,0,0))) -> myPlot
        
        ggsave(paste0("!01 Plots/01 Results/02 Num/Number_", gt, "_", a, ".png"), myPlot, 
               width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
        dev.off()
        
        tikz(file = paste0("!01 Plots/01 Results/02 Num/Number_", gt, "_", a, ".tikz"), 
             width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
        plot(myPlot)
        dev.off()
        
        # print(xtable(myResults,
        #              caption = "Power and FDR of each procedure",
        #              auto = TRUE),
        #       booktabs = TRUE)
    }
}
