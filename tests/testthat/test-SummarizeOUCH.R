## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: SummarizeOUCH")

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

### Define the SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=matrix(c(1,-1,0.5),ncol=1,nrow=3),
Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

### Now simulate the data.
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,NULL,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

## Here we do not do any recovery step
OUOU.summary<-SummarizeOUCH(phyltree,OUOUdata,OUOUparameters,regimes=NULL,t=c(1),dof=12)

pcmbase_OUOU_model_box<-PCMBase::PCM(model="OU",k=3)
pcmbase_OUOU_model_box$X0[] <- c(1,-1,0.5)
pcmbase_OUOU_model_box$Sigma_x[,, 1] <- OUOUparameters$Syy
pcmbase_OUOU_model_box$H[,, 1] <- OUOUparameters$A
pcmbase_OUOU_model_box$Theta[,1] <- c(1,-1,0.5)

testthat::expect_equivalent(PCMBase::PCMLik(t(OUOUdata), phyltree, pcmbase_OUOU_model_box, log = TRUE),OUOU.summary$t_1$LogLik)
## should be -7.171744349802688489603


