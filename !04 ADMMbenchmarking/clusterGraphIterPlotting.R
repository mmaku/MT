# Written by Micha³ Makowski

library(dplyr)
library(tidyr)
library(ggplot2)

load("./!02 Data/01 Binded/AllOneSimulations_Iter@18_10_18_19_58#88.RData")

a <- unique(finalResults$alpha)

plotFilename <- sprintf("./!01 Plots/clusterGraph/gslope_holm_glasso_block_diagonal_Iter_alpha_%s.png", a)

# names(finalResults)

finalResults %>% 
    dplyr::select(-c(penalizeDiagonal, localFDR)) %>%
    gather(FDR:Power, key = "metric", value = "value") %>% 
    group_by(n, p, alpha, graph.prob, graph.v, graph.u, procedure, metric) %>%
    ggplot(aes(x = iterations, y = value, color = metric)) +
    geom_point() +
    facet_grid(cols = vars(procedure), rows = vars(graph.u), labeller = label_context) 

    guides(color = guide_legend(override.aes = aes(size = 5)),
           shape = guide_legend(override.aes = aes(size = 5))) +
    labs(title = "Power and FDR for block diagonal matrices",
         subtitle = bquote(paste(alpha, " = ", .(a), ". Block size is 10.")),
         x = 'n') +
    theme(axis.title = element_text(angle = 0, inherit.blank = FALSE),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

ggsave(plotFilename, height = 5.5, width = 6)
