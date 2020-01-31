library(ggplot2)
library(dplyr)
#library(ggthemes)


df <- read.csv("Result/data_speedup.csv", colClasses = character(), sep = ",")

min <- df$t[1]
sp <- min / df$t

color <- as.character(df$v)

plot <- ggplot(df, aes(x = n, y = sp, col = color)) +
        labs(title = "", x = "number of process", y = "SpeedUp", color = color) +
        geom_line() + geom_point()

ggsave("pdf/speedup.pdf", plot=plot)
