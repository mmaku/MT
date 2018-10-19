# Written by Micha³ Makowski

library(dplyr)
library(tidyr)
library(ggplot2)

load("./!02 Data/01 Binded/AllOneSimulations_Prob@18_10_17_23_33#288.RData")

alphas <- unique(finalResults$alpha)
a <- alphas[1]
for(a in alphas)
{
    plotFilename <- sprintf("./!01 Plots/clusterGraph/gslope_holm_glasso_block_diagonal_probability_alpha_%s.png", alpha)
    finalResults %>% 
        filter(alpha == a) %>%
        dplyr::select(-c(penalizeDiagonal, iterations, SP, graph.u, graph.g, graphType)) %>%
        gather(FDR:Power, key = "metric", value = "value") %>% 
        group_by(n, p, alpha, graph.prob, procedure, metric) %>%
        ggplot(aes(x = factor(n), y = value, color = procedure, shape = metric)) +
        geom_jitter(width = 0.2, height = 0) +
        geom_hline(aes(yintercept = alpha)) +
        facet_wrap(vars(graph.v, graph.prob), labeller = label_context, scales = "fixed") +
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
}
