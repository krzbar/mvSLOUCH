## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: parametric.bootstrap")

library(mvSLOUCH)
library(PCMBase)


RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape. 
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(5)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),
Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))

### Now simulate the data.
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
BMdata<-BMdata[phyltree$tip.label,,drop=FALSE]

### Recover the parameters of the Brownian motion.
BMestim<-BrownianMotionModel(phyltree,BMdata)

### And finally obtain bootstrap confidence intervals for some parameters
BMbootstrap<-parametric.bootstrap(estimated.model=BMestim,phyltree=phyltree,
values.to.bootstrap=c("vX0","StS"),M.error=NULL,numboot=2)


testthat::expect_identical(length(BMbootstrap),2L)
testthat::expect_identical(length(BMbootstrap[[1]]),2L)
testthat::expect_identical(length(BMbootstrap[[1]][[1]]$bootstrapped.values),2L)
testthat::expect_identical(length(BMbootstrap[[1]][[2]]$bootstrapped.values),2L)

pcmbase_BM_model_box<-PCMBase::PCM(model="BM",k=3)
pcmbase_BM_model_box$X0[] <- BMbootstrap[[1]][[1]]$estimation.results$ParamsInModel$vX0[,1]
pcmbase_BM_model_box$Sigma_x[,, 1] <- BMbootstrap[[1]][[1]]$estimation.results$ParamsInModel$Sxx
testthat::expect_equivalent(PCMBase::PCMLik(t(BMbootstrap[[1]][[1]]$estimation.results$data), phyltree, pcmbase_BM_model_box, log = TRUE), BMbootstrap[[1]][[1]]$estimation.results$ParamSummary$LogLik)
## should be -5.963355179860567290007

testthat::expect_true(is.infinite(BMbootstrap[[1]][[2]]$estimation.results$ParamSummary$LogLik))
## should be -Inf


