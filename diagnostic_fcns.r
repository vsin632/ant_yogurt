diagnostics.plot<-function(mod.res){
  old.par = par(no.readonly = TRUE)
  par(mfrow=c(2, 2))
  par(mar=c(3, 3, 1, 0.5))
  hist(residuals(mod.res), probability=T, xlab="", ylab="", main="")
  mtext(text="histogram of residuals", side=3, line=0)
  x=seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out=100)
  lines(x, dnorm(x, mean=0, sd=sd(residuals(mod.res))))
  qqnorm(residuals(mod.res), main="", pch=19)
  qqline(residuals(mod.res))
  mtext(text="qq-plot of residuals", side=3, line=0)
  plot(fitted(mod.res), residuals(mod.res), pch=19)
  abline(h=0, lty=2)
  mtext(text="residuals against fitted values", side=3, line=0)
  par(old.par)
}

lev.thresh<-function(model.res){
	k=length(coefficients(model.res))
	n=length(residuals(model.res))
 return(2*(k+1)/n)
}

overdisp.test<-function(x){
  pr=residuals(x, type ="pearson")
  sum.dp=sum(pr^2)
  if(class(x)[[1]]=="mer"){
    if(length(grep(x=x@call, pattern="poisson"))==0){
      stop("oops, this isn't a model with poisson family... this function doesn't make sense")
    }
    xdf=length(residuals(x))-length(fixef(x))
  }else{
    if (x$family[[1]]!="poisson"){
      stop("oops, this isn't a model with poisson family... this function doesn't make sense")
    }
    xdf=length(residuals(x))-length(x$coefficients)
  }
  return(data.frame(chisq=sum.dp, df=xdf, P=1-pchisq(sum.dp, xdf), dispersion.parameter=sum.dp/xdf))
}
