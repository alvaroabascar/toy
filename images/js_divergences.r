library(ggplot2)
library(reshape)
a <- read.table('js_divergences.txt', header=TRUE)

png('./divergences.png', height=4, width=6, units='in', res=300)

a$pos = c(1:length(a$js))

p <- ggplot(data=a, aes(x=pos, y=js, fill=js)) +
     geom_bar(stat="identity") +
     xlab('Position') +
     ylab('Jensen-Shannon Divergence') +
     scale_x_continuous(breaks=c(1:9))

#p <- ggplot(data=a, aes(x=pos, y=js)) +
#     geom_bar() +
#     ylab('Jensen-Shannon divergence') +
#     xlab('position')

print(p)

dev.off()
