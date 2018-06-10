n = 20
B = 1000

gammahat = betahat = c()

for (i in 1:B) {
  y <- rweibull(n, 3, 100)
  gammahat[i]<-uniroot(function(x) n/x+sum(log(y))-n*
                      sum(y^x*log(y))/sum(y^x), c(1e-5,15))$root
  betahat[i]<- mean(y^gammahat)^(1/gammahat)
}

par(mfrow = c(1,2))
boxplot(gammahat, main = "gamma")
boxplot(betahat, main = "beta")

hist(gammahat, breaks = 40)
hist(betahat, breaks = 40)

par(mfrow = c(1,1))
