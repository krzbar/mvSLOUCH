## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: mvslouchModel")

library(mvSLOUCH)
library(PCMBase)


RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape. 
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(3)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

OUBMparameters<-list(vY0=matrix(1,ncol=1,nrow=1),A=matrix(0.5,ncol=1,nrow=1),
B=matrix(2,ncol=1,nrow=1),mPsi=matrix(1,1,1),
Syy=matrix(2,ncol=1,nrow=1),vX0=matrix(0,ncol=1,nrow=1),Sxx=diag(2,1,1),
Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1))


### Now simulate the data.
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,regimes=NULL,NULL)
OUBMdata<-OUBMdata[phyltree$tip.label,,drop=FALSE]

OUBMestim<-mvslouchModel(phyltree,OUBMdata,1,regimes=NULL,Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA="Positive",maxiter=c(1,2,1))


pcmbase_OUBM_model_box<-PCMBase::PCM(model="OU",k=2)
pcmbase_OUBM_model_box$X0[] <- c(OUBMestim$FinalFound$ParamsInModel$vY0[1,1],OUBMestim$FinalFound$ParamsInModel$vX0[1,1])
pcmbase_OUBM_model_box$Sigma_x[,, 1] <- rbind(cbind(OUBMestim$FinalFound$ParamsInModel$Syy,OUBMestim$FinalFound$ParamsInModel$Syx),cbind(OUBMestim$FinalFound$ParamsInModel$Sxy,OUBMestim$FinalFound$ParamsInModel$Sxx))
pcmbase_OUBM_model_box$H[,, 1] <- rbind(cbind(OUBMestim$FinalFound$ParamsInModel$A,OUBMestim$FinalFound$ParamsInModel$B),cbind(matrix(0,ncol=1,nrow=1),matrix(0,ncol=1,nrow=1)))
pcmbase_OUBM_model_box$Theta[,1] <- c(OUBMestim$FinalFound$ParamsInModel$mPsi[1,1],0)

testthat::expect_equivalent(PCMBase::PCMLik(t(OUBMdata), phyltree, pcmbase_OUBM_model_box, log = TRUE),OUBMestim$FinalFound$ParamSummary$LogLik)
## should be  -19.24322237182726169635



