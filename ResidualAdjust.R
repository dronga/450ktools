require(foreach)
require(doMC)

ResidualAdjust <- function(data, model, cores=1) {

  registerDoMC(cores=cores)
  
  results <- foreach(i=1:nrow(data),.combine=rbind) %dopar% {
    formula <- as.formula(paste("data[",i,",] ~ ", model))
    if(i%%1000==0)print(i)
    test <-try(residuals(lm(formula)))
    if(!inherits(test, "try-error") | length(test)!=ncol(data)){
        test
      } else { 
        rep(NA, ncol(data))
      }
  }
  rownames(results) <- rownames(data)
  return(results)
}
