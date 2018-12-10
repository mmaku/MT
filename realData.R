require(GEOquery, quietly = TRUE)
require(tibble, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(tidyr, quietly = TRUE)


data <- GDS2MA(getGEO(filename = "!02 Data/03 Real/GDS2771.soft.gz"))

correlation <- cor(t(data$M))



data$notes$channel_count

summary((dataSimple)
dataSimple$targets$disease.state
as_tibble(dataSimple$M) %>%
    cor()

cor((dataSimple$M))

dataSimple$M[,1]
