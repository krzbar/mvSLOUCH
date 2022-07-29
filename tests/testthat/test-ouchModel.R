## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: ouchModel")

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

OUOUparameters<-list(vY0=matrix(c(1,-1),nrow=2,ncol=1),
A=rbind(c(9,0),c(0,5)),mPsi=matrix(c(1,-1),ncol=1,nrow=2),
Syy=rbind(c(1,0.25),c(0,1)))

### Now simulate the data.
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes=NULL,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

OUOUestim<-ouchModel(phyltree,OUOUdata,regimes=NULL,Atype="SingleValueDiagonal",
Syytype="SingleValueDiagonal",diagA="Positive",maxiter=c(1,1))


pcmbase_OUOU_model_box<-PCMBase::PCM(model="OU",k=2)
pcmbase_OUOU_model_box$X0[] <- OUOUestim$MaxLikFound$ParamsInModel$vY0[,1]
pcmbase_OUOU_model_box$Sigma_x[,, 1] <- OUOUestim$MaxLikFound$ParamsInModel$Syy
pcmbase_OUOU_model_box$H[,, 1] <- OUOUestim$MaxLikFound$ParamsInModel$A
pcmbase_OUOU_model_box$Theta[,1] <- OUOUestim$MaxLikFound$ParamsInModel$mPsi[,1]

testthat::expect_equivalent(PCMBase::PCMLik(t(OUOUdata), phyltree, pcmbase_OUOU_model_box, log = TRUE),OUOUestim$MaxLikFound$LogLik)
## should be -240.4298491949068647955
testthat::expect_gt(OUOUestim$MaxLikFound$LogLik,OUOUestim$FinalFound$LogLik)
## should be -240.4298491949068647955 > -240.4298491949086837849
