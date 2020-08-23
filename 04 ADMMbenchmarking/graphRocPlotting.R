# Written by Micha³ Makowski

require(ggplot2, quietly = TRUE)
require(tikzDevice, quietly = TRUE)
require(xtable, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(tidyr, quietly = TRUE)
require(tibble, quietly = TRUE)

load("./!02 Data/01 Binded/01 AllOne/ROCfinal.RData")

result <- resultsList[[1]]
for(result in resultsList)
{
    SNR <- round((result$graphParameters$v+result$graphParameters$u)/result$graphParameters$v, 
                 digits = 2)
    gt <- result$graphType
    prob <- result$graphParameters$prob
    
    result$results %>%
        matrix(ncol = 3, byrow = TRUE) %>%
        as_tibble() %>%
        add_column(V4 = rep(dimnames(result$results)[[3]], 
                            each = dim(result$results)[2])) %>%
        mutate(method = recode_factor(V4,
                                      `banerjee.gLASSO` = "gLasso (Banerjee)", 
                                      `BH.gSLOPE` = "gSLOPE (BH)",
                                      `holm.gSLOPE` = "gSLOPE (Holm)")) %>%
        select(-V4) %>%
        group_by(method) -> result$finalResults
    
    colnames(result$finalResults) <- c(dimnames(result$results)[[1]], "method")
    
    result$finalResults %>%
        # select(-multiplier) %>%
        ggplot(aes(x = 1-SP, y = SN, color = method)) +
        geom_line(size = 1) +
        ylim(c(0,1)) +
        xlim(c(0,1)) +
        labs(subtitle = paste0("Setup: scaled $\\alpha$ = $0.05$, $100",
                               ifelse(gt != "scale-free", 
                                      paste0("$ variables, $", result$graphParameters$g, 
                                             "$ components, "), 
                                      "$ variables, "),
                               gt, " graph, $SNR$ = $", SNR,
                               ifelse(gt == "cluster", 
                                      paste0("$, $prob(x_{ij}\\neq0)=", 
                                             prob, "$."), 
                                      "$.")),
             y = "TPR",
             x = "FPR") +            
        scale_color_discrete(name = "Procedure:") +
        theme_bw(base_size = 8) +
        guides(color = guide_legend(override.aes = aes(size = 2)),
               shape = guide_legend(override.aes = aes(size = 2))) +
        theme(aspect.ratio = 1, 
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.box = "vertical",
              legend.spacing = unit(0, "cm"),
              plot.margin = margin(c(5,0,0,0)),
              legend.margin = margin(c(0,0,0,0))) -> myPlot
    
    ggsave(paste0("!01 Plots/01 Results/04 ROC/ROC_", gt, 
                  ifelse(gt == "cluster", paste0("_", prob, "_"), "_"), 
                  SNR, ".png"), myPlot, 
           width = 5.4, height = 5.4*myPlot$theme$aspect.ratio )
    dev.off()
    
    tikzTitle <- paste0("!01 Plots/01 Results/04 ROC/ROC_", gt, 
                        ifelse(gt == "cluster", paste0("_", prob, "_"), "_"), 
                        SNR, ".tikz")
    
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
    
    
}
