# Written by Micha³ Makowski

require(dplyr, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(tikzDevice, quietly = TRUE)
require(xtable, quietly = TRUE)


load("./!02 Data/01 Binded/01 AllOne/AllOneSimulations_HCSF_scaled_FDR_Corr@18_10_28_18_39#576.RData")

p    <- unique(bindedResults$p)
prob <- unique(bindedResults$graph.prob)
g    <- unique(bindedResults$graph.g)
gt   <- unique(bindedResults$graphType)[1]
a    <- unique(bindedResults$alpha)[1]


bindedResults %>%
    mutate(SNR = round((graph.v+graph.u)/graph.v, digits = 2)) %>%
    select(-c(SP, algIter, p, penalizeDiagonal, scaled, iterations, graph.g, graph.prob)) %>%
    as_tibble() -> bindedResults


bindedResults$SNR <- as.factor(bindedResults$SNR) 

labels1 <- sapply(levels(bindedResults$SNR), 
                  function(x) paste0("$SNR=", x,"$"))

# labels2 <- sapply(unique(bindedResults$graph.v),
#                  function(x) paste0("$\\textnormal{Off-diagonal}=", x,"$"))
# names(labels2) <- unique(bindedResults$graph.v)

for(gt in unique(bindedResults$graphType))
{
    for(a in unique(bindedResults$alpha))
    {
        
        bindedResults %>%
            filter(graphType == gt & alpha == a) %>%
            # filter(graph.v > 0) %>%
            mutate(procedure = recode_factor(procedure,
                                             `banerjee.gLASSO` = "gLasso (Banerjee)", 
                                             `BH.gSLOPE` = "gSLOPE (BH)",
                                             `holm.gSLOPE` = "gSLOPE (Holm)")) -> myResults
        
        gather(myResults, FDR:Power, key = "metric", value = "value") %>%
            ggplot(aes(x = factor(n), y = value, color = procedure, shape = metric)) +
            geom_jitter(width = 0.2, height = 0, size = 1.5) +
            geom_hline(aes(yintercept = a)) +
            ylim(c(0,1)) +
            # facet_wrap(vars(graph.v), labeller = labeller(graph.v = labels), 
            #            scales = "fixed", nrow = 4) +
            facet_wrap(vars(SNR), labeller = labeller(SNR = labels1), 
                       scales = "fixed", nrow = 2) +
            labs(subtitle = paste0("Setup: $\\alpha$ = $", a, "$, $",  p,
                                   ifelse(gt != "scale-free", 
                                          paste0("$ variables, $", g, "$ components, "), 
                                          "$ variables, "),
                                   gt, " graph",
                                   ifelse(gt == "cluster", 
                                          paste0(", $prob(x_{ij}\\neq0)=", prob, "$."), 
                                          ".")),
                 y = "Value",
                 x = "Sample size $n$") +
            scale_color_discrete(name = "Procedure:") +
            scale_shape_discrete(name = "Measure:") +
            theme_bw(base_size = 8) +
            guides(color = guide_legend(override.aes = aes(size = 2)),
                   shape = guide_legend(override.aes = aes(size = 2))) +
            # theme(aspect.ratio = 8/16, 
            #       plot.margin = margin(c(0,0,0,0)),
            #       legend.margin = margin(c(0,10,0,0))) -> myPlot
            theme(aspect.ratio = 12/16, 
                  legend.position = "bottom",
                  legend.direction = "horizontal",
                  legend.box = "vertical",
                  legend.spacing = unit(0, "cm"),
                  plot.margin = margin(c(5,0,0,0)),
                  legend.margin = margin(c(0,0,0,0))) -> myPlot
        
        ggsave(paste0("!01 Plots/01 Results/03 Corr/Corr_", gt, "_", a, ".png"), myPlot, 
               width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
        dev.off()
        
        tikzTitle <- paste0("!01 Plots/01 Results/03 Corr/Corr_", gt, "_", a, ".tikz") 
        
        tikz(file = tikzTitle, 
             width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
        plot(myPlot)
        dev.off()
        
        lines <- readLines(con = tikzTitle)
        lines <- lines[-which(grepl("\\path\\[clip\\]*", x = lines, perl=F))]
        lines <- lines[-which(grepl("\\path\\[use as bounding box*", x = lines, perl=F))]
        lines <- gsub(pattern = "SNR", replace = "\\SNR", x = lines, fixed = TRUE)
        lines <- gsub(pattern = "prob", replace = "\\prob", x = lines, fixed = TRUE)
        writeLines(lines,con = tikzTitle)
        
        # print(xtable(myResults,
        #              caption = "Power and FDR of each procedure",
        #              auto = TRUE),
        #       booktabs = TRUE)
    }
}
