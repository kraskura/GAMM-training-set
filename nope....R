# trash BELOW

#### 3. Sole egg data



```{r,echo=FALSE}

data(sole)
data(coast)
sole$off <- log(sole$a.1-sole$a.0)# model offset term
sole$a<-(sole$a.1+sole$a.0)/2 # mean stage age
solr<-sole # make copy for rescaling
solr$station <- factor(with(solr, paste (-la, -lo, -t, sep="")))

{
  par(mfrow=c(2,3))
  sample.t <- unique(sole$t)
  stage <- 1
  for (i in 1:5){
    egg<-sole[sole$stage==stage & sole$t==sample.t[i],]
    plot(egg$lo,egg$la,xlab="lo",ylab="la",main=paste("day",sample.t[i]), cex=egg$eggs/4, xlim=range(sole$lo), ylim=range(sole$la), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
    points(egg$lo,egg$la, pch=".", col=2)
    lines(coast)
    # points(solr[solr$stage==solr & solr$t==sample.t[i],], solr$lo,solr$la,col="blue")
  }
  
  ## boundary definition list and knots suitable for soap film smoothing
  bnd <- list(list(lo=c(-6.74,-5.72,-5.7 ,-5.52,-5.37,-5.21,-5.09,-5.02,
  -4.92,-4.76,-4.64,-4.56,-4.53,-4.3,-4.16,-3.8 ,-3.8,-5.04,-6.76,
  -6.74),
  la=c(50.01,50.02,50.13,50.21,50.24,50.32,50.41,50.54,50.59,50.64,
  50.74,50.86,51.01,51 ,51.2,51.22,51.61,51.7,51.7,50.01)))
  knt <- list(lo=c(-4.643,-5.172,-5.638,-6.159,-6.665,-6.158,-5.656,-5.149,
  -4.652,-4.154,-3.901,-4.146,-4.381,-4.9,-5.149,-5.37,-5.866,-6.36,-6.635,
  -6.12,-5.626,-5.117,-4.622,-4.695,-4.875,-5.102,-5.609,-5.652,-5.141,
  -5.354,-5.843,-6.35,-6.628,-6.127,-5.63,-5.154,-5.356,-5.652,-5.853,
  -6.123),
  la=c(51.626,51.61,51.639,51.638,51.376,51.377,51.373,51.374,51.374,
  51.376,51.379,51.226,51.129,51.194,51.083,51.147,51.129,51.151,50.901,
  50.891,50.959,50.958,50.942,50.728,50.676,50.818,50.825,50.684,50.693,
  50.568,50.564,50.626,50.397,50.451,50.443,50.457,50.325,50.193,50.322,
  50.177))
  points(knt$lo,knt$la,pch=19,col=2,cex=.6)
  lines(bnd[[1]]$lo,bnd[[1]]$la,col=2)
  
  
}

solr.dummy<-expand.grid(t=c(50,68, 86, 104, 122, 140), lo=c(seq(-6.5, -4, 0.5)), la=c(seq(50, 51.7, 0.5)), a=0, off=0)
# solr.dummy$station <- factor(with(solr.dummy, paste (-la, -lo, -t, sep="")))

# solr.dummy$pred.b4<-predict(b4,newdata = solr.dummy)
# solr$pred.b4<-predict(b4)
  
```



##### 3.1 Generalized model



##### 3.2 GAM 

```{r mgcv, eval=FALSE, include=FALSE}

som<-gamm(eggs~te(lo, la, t, bs=c("tp", "tp"), k=c(25,5), d=c(2,1)) + s(t, k=5, by=a) + offset(off), family=quasipoisson, data=solr, random=list(station=~1))

solr.dummy$pred.som<-predict(som$gam,newdata = solr.dummy)

```


##### 3.3 BAM?! yes :) 

Faster alternative: 
(this used REML for working model estimation

```{r mgcv, echo=FALSE, include=FALSE}
som1<-bam(eggs~te(lo, la, t, bs=c("tp", "tp"), k=c(25,5), d=c(2,1)) + s(t, k=5, by=a) + offset(off) + s(station, bs="re"), family=quasipoisson, data=solr)

solr.dummy$som1<-predict(som1, exclude="s(station)", newdata = solr.dummy)


```


``` {r echo=FALSE}
  {
  par(mfrow=c(2,3))
  sample.t <- unique(solr.dummy$t)
  # sample.t <- unique(sole$t)
  for (i in 1:6){
    egg.pred<-solr.dummy[solr.dummy$t==sample.t[i],]
    plot(egg.pred$lo, egg.pred$la,xlab="lo",ylab="la",main=paste("day",sample.t[i]), cex=egg.pred$pred.som/4, xlim=range(sole$lo), ylim=range(sole$la), cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
    points(egg.pred$lo,egg.pred$la, pch=".", col=2)
    lines(coast)
    # lines(solr.dummy[solr.dummy$pred.som== & solr$t==sample.t[i],], solr$lo,solr$la,col="blue", add=TRUE)
  }
  
 
  # points(knt$lo,knt$la,pch=19,col=2,cex=.6)
  # lines(bnd[[1]]$lo,bnd[[1]]$la,col=2)
  
  
}
```


```{r gamair, echo=FALSE }

## 3.3.5
data(sole)
sole$off <- log(sole$a.1-sole$a.0)# model offset term
sole$a<-(sole$a.1+sole$a.0)/2 # mean stage age
solr<-sole # make copy for rescaling
solr$t<-solr$t-mean(sole$t)
solr$t<-solr$t/var(sole$t)^0.5
solr$la<-solr$la-mean(sole$la)
solr$lo<-solr$lo-mean(sole$lo)
b <- glm(eggs ~ offset(off)+lo+la+t+I(lo*la)+I(lo^2)+I(la^2)
+I(t^2)+I(lo*t)+I(la*t)+I(lo^3)+I(la^3)+I(t^3)+
I(lo*la*t)+I(lo^2*la)+I(lo*la^2)+I(lo^2*t)+
I(la^2*t)+I(la*t^2)+I(lo*t^2)+ a +I(a*t)+I(t^2*a),
family=quasi(link=log,variance="mu"),data=solr)
summary(b)
b1 <- update(b, ~ . - I(lo*t))
b4 <- update(b1, ~ . - I(lo*la*t) - I(lo*t^2) - I(lo^2*t))
anova(b,b4,test="F")
par(mfrow=c(1,2)) # split graph window into 2 panels
plot(fitted(b4)^0.5,solr$eggs^0.5) # fitted vs. data plot
plot(fitted(b4)^0.5,residuals(b4)) # resids vs. sqrt(fitted)

plot(predict(b4)~solr$egg)

solr$perdictedb4<-predict(b4)

## 3.5.1
rf <- residuals(b4,type="d") # extract deviance residuals
## create an identifier for each sampling station
solr$station <- factor(with(solr,paste(-la,-lo,-t,sep="")))
## is there evidence of a station effect in the residuals?
solr$rf <-rf
rm <- lme(rf~1,solr,random=~1|station)
rm0 <- lm(rf~1,solr)
anova(rm,rm0)
## following is slow...




```