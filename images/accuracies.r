library(ggplot2)

a <- read.table('accuracies.txt', header=TRUE)

a$model <- factor(a$model, levels=c("toy", "u1", "u1_all", "tia1", "tia2", "tia3"))

png('./accuracies.png', height=4, width=6, units='in', res=300)

p <- ggplot(data=a, aes(x=model, y=sensitivity, fill=sensitivity)) +
     geom_bar(stat="identity") +
     xlab('Model') +
     ylab('Sensitivity (%)')

print(p)

dev.off()
