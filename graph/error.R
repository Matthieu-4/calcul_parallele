library(ggplot2)
library(dplyr)
#library(ggthemes)


df <- read.csv("Result/data_speedup.csv", colClasses = character(), sep = ",")


version <- as.character(df$v)
x1 <- df$n[1:3] - 0.25
x2 <- df$n[4:6]
x3 <- df$n[7:9] + 0.25
x_off <- c(x1, x2, x3)



plot <- ggplot(df, aes(x = x_off, y = e, col = version)) +
        labs(title = "", x = "values of h", y = "Numerical error") +
        geom_bar(stat="identity", aes(fill=version))

ggsave("pdf/error.pdf", plot=plot)
