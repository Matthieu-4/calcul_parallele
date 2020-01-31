library(ggplot2)
library(dplyr)
#library(ggthemes)


df <- read.csv("Result/data_h.csv", colClasses = character(), sep = ",")


version <- as.character(df$v)
x1 <- df$h[1:10] - 0.25
x2 <- df$h[11:20]
x3 <- df$h[21:30] + 0.25
x_off <- c(x1, x2, x3)
print(df$h[1:6])
print(x1)
print(x_off)
plot <- ggplot(df, aes(x = x_off, y = t, col = version)) +
        labs(title = "", x = "values of h", y = "Time (us)") +
        #geom_line() + geom_point() +
        geom_bar(stat="identity", aes(fill=version))

ggsave("pdf/h.pdf", plot=plot)


plot <- ggplot(df, aes(x = x_off, y = e, col = version)) + ylim(0, 2.5e-157)+
        labs(title = "", x = "values of h", y = "Numerical error") +
        #geom_line() + geom_point() +
        geom_bar(stat="identity", aes(fill=version))

ggsave("pdf/h_error.pdf", plot=plot)
