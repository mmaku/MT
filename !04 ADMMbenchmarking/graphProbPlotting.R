# Written by Micha³ Makowski

require(dplyr, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(tikzDevice, quietly = TRUE)
require(xtable, quietly = TRUE)


load("./!02 Data/01 Binded/01 AllOne/AllOneSimulations_Cl_scaled_FDR_Prob@18_10_28_00_36#288.RData")

p    <- unique(bindedResults$p)
prob <- unique(bindedResults$graph.prob)[1]
g    <- unique(bindedResults$graph.g)
gt   <- unique(bindedResults$graphType)
a    <- unique(bindedResults$alpha)[1]

bindedResults %>%
    mutate(SNR = round((graph.v+graph.u)/graph.v, digits = 2)) %>%
    select(-c(SP, algIter, p, graphType, penalizeDiagonal, scaled, iterations, graph.g)) %>%
    as_tibble() -> bindedResults

# bindedResults %>%
#     select(-c(localFDR, penalizeDiagonal, iterations, SP)) %>%
#     filter(procedure != "gLASSO") %>%
#     as_tibble() -> bindedResults

bindedResults$SNR <- as.factor(bindedResults$SNR) 

rowLabels <- sapply(levels(bindedResults$SNR), 
                  function(x) paste0("$SNR=", x,"$"))

colLabels <- sapply(unique(bindedResults$graph.prob),
                    function(x) paste0("$prob(x_{ij}\\neq0)=", x, "$"))
names(colLabels) <- unique(bindedResults$graph.prob)


for(a in unique(bindedResults$alpha))
{
    bindedResults %>%
        filter(alpha == a) %>%
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
        facet_grid(cols = vars(graph.prob), rows = vars(SNR), 
                   labeller = labeller(.rows =  rowLabels,
                                       .cols = colLabels), 
                   scales = "fixed") +
        labs(subtitle = paste0("Setup: $\\alpha$ = $", a, "$, $",  
                               p, "$ variables, $", g, "$ components, ",
                               gt, " graph."),
             y = "Value",
             x = "Sample size $n$") +
        scale_color_discrete(name = "Procedure:") +
        scale_shape_discrete(name = "Measure:") +
        theme_bw(base_size = 8) +
        guides(color = guide_legend(override.aes = aes(size = 2)),
               shape = guide_legend(override.aes = aes(size = 2))) +
        theme(aspect.ratio = 12/16, 
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.box = "vertical",
              legend.spacing = unit(0, "cm"),
              plot.margin = margin(c(5,0,0,0)),
              legend.margin = margin(c(0,0,0,0))) -> myPlot
    
    ggsave(paste0("!01 Plots/01 Results/01 Prob/Prob_", a, ".png"), myPlot, 
           width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
    dev.off()
    
    tikzTitle <- paste0("!01 Plots/01 Results/01 Prob/Prob_", a, ".tikz") 
    
    tikz(file = tikzTitle, 
         width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
    plot(myPlot)
    dev.off()
    
    lines <- readLines(con = tikzTitle)
    lines <- lines[-which(grepl("\\path\\[clip\\]*", x = lines, perl=F))]
    lines <- lines[-which(grepl("\\path\\[use as bounding box*", x = lines, perl=F))]
    lines <- gsub(pattern = "SNR", replace = "\\SNR", x = lines, fixed = TRUE)
    lines <- gsub(pattern = "prob", replace = "\\prob", x = lines, fixed = TRUE)
    writeLines(lines, con = tikzTitle)
    
    # print(xtable(myResults,
    #              caption = "Power and FDR of each procedure",
    #              auto = TRUE),
    #       booktabs = TRUE)
}
