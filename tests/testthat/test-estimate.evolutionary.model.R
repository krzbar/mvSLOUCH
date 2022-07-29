## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: estimate.evolutionary.model")

library(mvSLOUCH)
library(PCMBase)


RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape.
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(4)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

### Define a vector of regimes.
regimes<-c("small","small","large","small","large","small")

### Define SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),"large"=c(-1,1,0.5)),
Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

### Now simulate the data.
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

## set up for a trivial, single model setup case (for running time)
## in a real analysis you should carefully choose between what models
## you want to do model selection
model_setups<-list(list(evolmodel="bm"))

### Try to recover the parameters of the OUOU model.
### maxiter here set to minimal working possibility, in reality it should be larger
### e.g. default of c(10,50,100)
estimResults<-estimate.evolutionary.model(phyltree,OUOUdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=1,model.setups=model_setups,predictors=c(3),
kY=2,doPrint=TRUE,pESS=NULL,maxiter=c(1,1,1))


pcmbase_BM_model_box<-PCMBase::PCM(model="BM",k=3)
pcmbase_BM_model_box$X0[] <- estimResults[[1]]$BestModel$ParamsInModel$vX0[,1]
pcmbase_BM_model_box$Sigma_x[,, 1] <-  estimResults[[1]]$BestModel$ParamsInModel$Sxx
testthat::expect_equivalent(PCMBase::PCMLik(t(OUOUdata), phyltree, pcmbase_BM_model_box, log = TRUE),estimResults[[1]]$BestModel$LogLik)

## should be -13.57518237951901873828



