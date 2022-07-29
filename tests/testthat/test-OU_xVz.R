## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: OU_xVz")

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

### Define a vector of regimes.
regimes<-c("small","small","large","small","small","large","large","large")

### Define SDE parameters to be able to simulate data under the OUOU model.
## 3D model
## OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
## A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),"large"=c(-1,1,0.5)),
## Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))
## 2D model used to reduce running time on CRAN
OUOUparameters<-list(vY0=matrix(c(1,-1),nrow=2,ncol=1),
A=rbind(c(9,0),c(0,5)),mPsi=cbind("small"=c(1,-1),"large"=c(-1,1)),
Syy=rbind(c(1,0.25),c(0,1)))

### Now simulate the data.
OUOUdata1<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata1<-OUOUdata1[phyltree$tip.label,,drop=FALSE]

OUOUdata2<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata2<-OUOUdata2[phyltree$tip.label,,drop=FALSE]

xVz_output<-OU_xVz(OUOUdata1, OUOUdata2, phyltree, OUOUparameters, M.error=NULL, do_centre="evolutionary_model", regimes = regimes)

testthat::expect_equivalent(xVz_output$xVz, -90.31009696162655586704)



