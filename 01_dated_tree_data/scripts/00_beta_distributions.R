x <- 0:1000/1000

pdf('../outputs/priors.pdf', 6, 5)
plot(x, dbeta(x, 1.5, 2.5), type = 'l', lty = 1, col = 'steelblue', lwd = 3, ylab = 'f(x)')
lines(x, dbeta(x, 1.2, 1.2), type = 'l', lty = 4, col = 'darkolivegreen', lwd = 3)
title('Beta probability density functions')
dev.off()

library(ggplot2)
theme_set(theme_bw())

pdf('../outputs/priors_gg.pdf', 5.4/1.5, 3.75/1.5)
ggplot() +
  geom_line(aes(x, dbeta(x, 1.5, 2.5)), lty = 1, col = 'steelblue', lwd = 1) +
  geom_line(aes(x, dbeta(x, 1.2, 1.2)), lty = 4, col = 'darkolivegreen', lwd = 1) +
  ylab('Density')
dev.off()
