###############################################################################
### a function to simulate (0,1) genotype markers with tunable
### linkage disequilibrium
### n=number of individuals
### nmarkers=number of markers
### recombprob=probablities of recombination; at each successive locus, one of 
### these values is sampled randomly to determine how many of the individuals
### will "flip" their genotype; the more '0' values in this vector, the more
### LD there will be in the markers.  Higher values will result in more 
### recombination.
###############################################################################

simgeno = function(n,nmarkers=1000,recombprob=c(0,0,0,0,0.01,0.3)){
 geno = matrix(0,n,nmarkers)
 xn = sample(c(0,1),n,replace=T)
 geno[,1] = xn
 draw = ceiling(recombprob*n)
 for(i in 2:nmarkers){
  these = sample(1:n,sample(draw,1))
  xn[these] = !xn[these]
  geno[,i] = xn   
 }
 return(geno)
}

###############################################################################
### a function to simulate traits, given an input logical vector and an effect
### size for those "with" the trait (i.e. TRUE)
### logic=a logical vector indicating which samples "have" the trait (TRUE) and
### which do not (FALSE)
### effect=effect size (separation between those that have the trait and those
### that don't
### spread=value that controls the spread of the distributions 
###############################################################################

simtrait = function(logic,effect=1,spread=0.1){
 spread = spread*abs(effect)
 logic = as.logical(logic)
 out = numeric(length(logic))
 out[!logic] = rnorm(sum(!logic),0,spread)
 out[logic] = rnorm(sum(logic),effect,spread)
 return(out)
}


##############################################################################
### function for extracting selection frequencies from an RF
### rf=the RF from which selection frequencies are desired
##############################################################################
rfsf = function(rf){
 vu = randomForest::varUsed(rf)
 sf = vu/sum(vu)
 names(sf) = rownames(rf$importance)
 return(sf)
}

##############################################################################
### function to estimate the selection bias in RF when the null hypothesis is
### true (i.e. no association between response and predictors)
### x=predictor matrix (genotype matrix)
### ntree=number of trees 
### verbose=print progress?
### mult=create a multiplicative bias correction? (otherwise additive)
### ...=additional arguments to be passed to 'randomForest'
##############################################################################
estBias = function(x,ntree,verbose=TRUE,mult=FALSE,...){
 if(ntree%%10!=0) stop("ntree must be a multiple of 10")
 vu = numeric(ncol(x))
 ni = ntree/10
 for(i in 1:ni){
  rf = randomForest::randomForest(y=rnorm(nrow(x)),x=x,ntree=10,...)
  vu = vu + randomForest::varUsed(rf)
  if(i%%10==0 && verbose){
   cat(round(100*i/ni,0),"percent complete")
   cat("\n")
  }
 }
 vu = vu/sum(vu)
 if(mult){
  return(mean(vu)/vu)
 }else{
  return(vu-mean(vu))
 }
}

######################################################################
### plot density with data points
######################################################################
pdens = function(x,...){
 plot(density(x),...)
 points(x,rep(0,length(x)),col='red',pch="|")
}

?parApply
library(snow)

head(geno)

head(pheno)
pheno


