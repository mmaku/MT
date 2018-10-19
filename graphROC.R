
results <- graphROC(n = 150, 
                    p = 200, 
                    graphType = "cluster",
                    graphParameters = NULL, 
                    method = "banerjee",
                    alpha = 0.05, 
                    alphaMulti = c(1/30, 3/2),
                    alphaLength = 15,
                    penalizeDiagonal = FALSE, 
                    iterations = 5,  # Numer of graphs simulated to calculated FDR
                    epsilon = 10e-4, 
                    verbose = TRUE)


colfunc<-colorRampPalette(c("black","yellow"))
plot(1-results$SP, results$SN, type = "b", xlim = c(0,1), ylim = c(0,1), pch = seq_along(results$multi))
legend("bottomright", legend = "multi = 1", pch = which.min(abs(results$multi - 1)))


