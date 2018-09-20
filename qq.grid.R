qq.ci <- function(obs,main="",col="black") {

obs <- -log10(obs)
obs <- sort(obs)
N <- length(obs) ## number of p-values

## create the null distribution
## (-log10 of the uniform)
null <- -log(1:N/N,10)
MAX <- max(c(obs,null),na.rm=T)




## create the confidence intervals
c95 <- rep(0,N)
c05 <- rep(0,N)

## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury)

for(i in 1:N){
c95[i] <- qbeta(0.95,i,N-i+1)
c05[i] <- qbeta(0.05,i,N-i+1)
}

## plot the two confidence lines
par(mar=c(5,5,4,2))
#plot(null, -log(c95,10), ylim=c(0,max(obs,na.rm=T)),xlim=c(0,max(null,na.rm=T)), main=main, type="l",
#axes=T, col="grey", xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), cex.lab=1.3, cex.axis=1.3)
#par(new=T)
plot(null, -log(c05,10),  ylim=c(0,max(obs,na.rm=T)),xlim=c(0,max(null,na.rm=T)), main=main,type="l",
axes=T, col="grey", xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), cex.lab=1.3, cex.axis=1.3)

polygon(c(null,rev(null)),c(-log(c05,10),rev(-log(c95,10))),col="grey",border=NA)
## add the diagonal
abline(0,1,col="red",lwd=2)
#par(new=T)

grid(lwd=3)

## add the qqplot
points(sort(null),sort(obs), ylim=c(0,max(null,na.rm=T)),xlim=c(0,max(obs,na.rm=T)), xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),main=main, pch=19, cex=0.8, col=col)

}

qq.ci.add <- function(obs,main="",col="red") {
obs <- -log10(obs)
obs <- sort(obs)
N <- length(obs) ## number of p-values
null <- -log(1:N/N,10)
points(sort(null),sort(obs), ylim=c(0,max(null,na.rm=T)),xlim=c(0,max(obs,na.rm=T)), xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),main=main, pch=19, cex=0.8, col=col)
}
