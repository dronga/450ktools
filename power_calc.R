
## Power Calculations for EWAS

library(pwr)

n.meth.psets = 850000
alpha = .05
#case/con effect size detectable at N=500 
out.cascon = pwr.t.test(n = 500, sig.level = alpha / n.meth.psets, power = .8, type = "two.sample")
paste("Powered to detect standardized effect size of", round(out.cascon$d, 2))


#power curves for N=100, N=500, N=1000 

x <- seq(0.01,0.6,by=0.01)
y100 <- sapply(x,  function(z){pwr.t.test(n = 100, sig.level = alpha / n.meth.psets, d = z, type = "two.sample")$power})
y500 <- sapply(x, function(z){pwr.t.test(n = 500, sig.level = alpha / n.meth.psets, d = z, type = "two.sample")$power})
y1000 <- sapply(x, function(z){pwr.t.test(n = 1000, sig.level = alpha / n.meth.psets, d = z, type = "two.sample")$power})

par(mar=c(5,6,4,2))
plot(x,y1000,type="l", xlab="Standardized effect size", ylab="Power", main="Case/Control association", col="purple", cex.lab=1.3, cex.axis=1.3, lwd=3)
lines(x,y500,xlab="Standardized effect size", ylab="Power", main="Case/Control association", col="orange", cex.lab=1.3, cex.axis=1.3, lwd=3)
lines(x,y100,xlab="Standardized effect size", ylab="Power", main="Case/Control association, N = 500", col="cornflowerblue", cex.lab=1.3, cex.axis=1.3, lwd=3)
grid(lwd=3)

legend("topleft", legend=c("N =  100", "N =  500", "N = 1000"), col=c("cornflowerblue", "orange", "purple"), lty=1, lwd=3, bg="white")
