## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: SummarizeMVSLOUCH")

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

### Define SDE parameters to be able to simulate data under the mvOUBM model.
OUBMparameters<-list(vY0=matrix(1,ncol=1,nrow=1),A=matrix(0.5,ncol=1,nrow=1),
B=matrix(c(2),ncol=1,nrow=1),mPsi=matrix(1,ncol=1,nrow=1),
Syy=matrix(2,ncol=1,nrow=1),vX0=matrix(0,ncol=1,nrow=1),Sxx=diag(1,1,1),
Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1))

### Now simulate the data.
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,NULL,NULL)
OUBMdata<-OUBMdata[phyltree$tip.label,,drop=FALSE]

## Here we do not do any recovery step
OUBM.summary<-SummarizeMVSLOUCH(phyltree,OUBMdata,OUBMparameters,regimes=NULL,t=c(1),dof=6,maxiter=2)

pcmbase_OUBM_model_box<-PCMBase::PCM(model="OU",k=2)
pcmbase_OUBM_model_box$X0[] <- c(1,0)
pcmbase_OUBM_model_box$Sigma_x[,, 1] <- rbind(cbind(OUBMparameters$Syy,OUBMparameters$Syx),cbind(OUBMparameters$Sxy,OUBMparameters$Sxx))
pcmbase_OUBM_model_box$H[,, 1] <- rbind(cbind(OUBMparameters$A,OUBMparameters$B),cbind(matrix(0,ncol=1,nrow=1),matrix(0,ncol=1,nrow=1)))
pcmbase_OUBM_model_box$Theta[,1] <- c(1,0)

testthat::expect_equivalent(PCMBase::PCMLik(t(OUBMdata), phyltree, pcmbase_OUBM_model_box, log = TRUE),OUBM.summary$t_1$LogLik)
## should be  -10.4025115800806826627



