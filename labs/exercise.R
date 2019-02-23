# andare a vedere la lezione di Egidi

library(faraway)
library(mgcv)
colnames(ozone)

ozone.fit <- gam(log(O3) ~s(vh)+s(wind)+s(humidity)+s(temp)+s(ibh)+s(dpg)+s(ibt)+s(vis)+s(doy), 
                 data=ozone, family=gaussian)
summary(ozone.fit)

par(mfrow=c(3,3))
plot(ozone.fit)
par(mfrow=c(1,1))

ozone.gam <- gam(log(O3) ~vh+wind+s(humidity)+s(temp)+s(ibh)+s(dpg)+ibt+s(vis)+s(doy), 
                 data=ozone, family=gaussian)

summary(ozone.gam)


ozone.gam.gamma <- gam(O3 ~vh+wind+s(humidity)+s(temp)+s(ibh)+s(dpg)+ibt+s(vis)+s(doy), 
                  data=ozone, family=Gamma(link=log))

ozone.glm <- glm(O3 ~vh+wind+humidity+temp+ibh+dpg+ibt+vis+doy, 
                  data=ozone, family=Gamma(link=log))

AIC(ozone.gam, ozone.gam.gamma, ozone.glm)
# il modello migliore è il gam, quello con AIC più basso


# 2 - PROVIDE THE OPTIMAL DEGREE OF SMOOTHING
sp <- ozone.gam$sp

tuning.scale<-c(1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,
                1e3,1e4,1e5) 
scale.exponent <- log10(tuning.scale) 
n.tuning <- length(tuning.scale) 
edf <- rep(NA,n.tuning)  
min2ll <- rep(NA,n.tuning) 
aic <- rep(NA,n.tuning)  


for (i in 1:n.tuning) {   
  gamobj <- gam(Volume ~ s(Height) + s(Girth), family=Gamma(link=log),
                data=trees, sp=tuning.scale[i]*sp) 
  min2ll[i]<--2*logLik(gamobj) 
  edf[i]<-sum(gamobj$edf)+1  
  aic[i]<-AIC(gamobj)
}

par(mfrow=c(2,2)) 
plot(scale.exponent,min2ll,type="b",main="-2 log likelihood") 
plot(scale.exponent,edf,ylim=c(0,70),type="b",main="effective number of parameters")  
plot(scale.exponent,aic,type="b",main="AIC") 

# find the minimum aic
min.aic<-1e100 
opt.tuning.scale<-NULL 
for (i in 1:n.tuning) { 
  if (aic[i]<min.aic) { 
    min.aic<-aic[i] 
    opt.tuning.scale<-tuning.scale[i] 
  } 
} 
opt.sp<-opt.tuning.scale*sp

# find the minimum bic

# fitting the final model with the optimal level of smoothing

gam.2.opt <- gam(Volume ~ s(Height)+s(Girth), family=Gamma(link=log),data=trees, sp=opt.sp)

AIC(gam.2, gam.2.opt)