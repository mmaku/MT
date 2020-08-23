# Written by Micha³ Makowski

require(dplyr, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(tikzDevice, quietly = TRUE)
require(xtable, quietly = TRUE)


load("./!02 Data/01 Binded/01 AllOne/AllOneSimulations_HC_scaled_FDR_Num@18_11_01_16_40#288.RData")

finalResults %>%
    select(-c(penalizeDiagonal, iterations, SP)) %>%
    as_tibble() -> finalResults1

load("./!02 Data/01 Binded/01 AllOne/AllOneSimulations_Cl_scaled_FDR_Prob@18_10_28_00_36#288.RData")

bindedResults %>%
    select(-c(penalizeDiagonal, iterations, SP)) %>%
    filter(graph.prob == 0.5) %>%
    as_tibble() -> finalResults2

load("./!02 Data/01 Binded/01 AllOne/AllOneSimulations_HCSF_scaled_FDR_Corr@18_10_28_18_39#576.RData")

bindedResults %>%
    select(-c(penalizeDiagonal, iterations, SP)) %>%
    filter(graphType == "hub" & ((graph.v == 0.7)|(graph.v == 0.3)|graph.v == -1)) %>%
    as_tibble() -> finalResults3


finalResults <- rbind(finalResults1, finalResults2, finalResults3)

p    <- unique(finalResults$p)
prob <- unique(finalResults$graph.prob)
g    <- unique(finalResults$graph.g)[1]
gt   <- unique(finalResults$graphType)[1]
a    <- unique(finalResults$alpha)[1]

finalResults %>%
    mutate(SNR = round((graph.v+graph.u)/graph.v, digits = 2)) %>%
    select(-c(algIter, p, scaled, graph.v, graph.u, graph.prob)) -> finalResults

finalResults$SNR <- as.factor(finalResults$SNR) 

rowLabels <- sapply(levels(finalResults$SNR), 
                    function(x) paste0("$SNR=", x,"$"))

colLabels <- sapply(unique(finalResults$graph.g),
                    function(x) paste0("$\\#[\\textnormal{conectivity components}] = ", x, "$"))
names(colLabels) <- unique(finalResults$graph.g)

for(gt in unique(finalResults$graphType))
{
    for(a in unique(finalResults$alpha))
    {
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
            facet_grid(cols = vars(graph.g), rows = vars(SNR), 
                       labeller = labeller(.rows =  rowLabels,
                                           .cols = colLabels), 
                       scales = "fixed") +
            labs(subtitle = paste0("Setup: $\\alpha$ = $", a, "$, $",  p,
                                   "$ variables, ", gt, " graph",
                                   ifelse(gt == "cluster", 
                                          paste0(", $prob(x_{ij}\\neq0)=", prob, "$."), 
                                          ".")),
                 y = "Value",
                 x = "Sample size $n$") +
            scale_color_discrete(name = "Procedure:") +
            scale_shape_discrete(name = "Measure:") +
            # theme_bw() +
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
        
        ggsave(paste0("!01 Plots/01 Results/02 Num/Number_", gt, "_", a, ".png"), myPlot, 
               width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
        dev.off()
        
        tikzTitle <- paste0("!01 Plots/01 Results/02 Num/Number_", gt, "_", a, ".tikz")
        
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
