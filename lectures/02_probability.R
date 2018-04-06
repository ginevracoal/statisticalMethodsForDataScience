library(MASS)

# BINOMIAL DISTRIBUTION
x <- rbinom(10^5, size = 1, prob = 0.7)

# NORMAL DISTRIBUTION
# rnorm generates random deviates
rnorm(10) # If mean or sd are not specified they assume the 
          # default values of 0 and 1
rnorm(10, m=5) # m = mean
rnorm(10, m=5, s=2) # s = standard deviation
y <- rnorm(10^5, m = x * 5, s = 1) ### Y| X = x ~ N(x * 5, 1)

# density function graph on the interval [-4,10]
hist.scott(y, main = "", xlim = c(-4, 10))

xx <- seq(-4, 10, l = 1000) # l is the length of the sequence

# dnorm gives the density
ff <- 0.3 * dnorm(xx, 0) + 0.7 * dnorm(xx, 5)

### This is a mixture of normal distributions
hist.scott(y, main = "", xlim = c(-4, 10))
lines(xx, ff, col = "red", lwd = 2)
