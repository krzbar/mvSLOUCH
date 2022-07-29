## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: fitch.mvsl")

library(mvSLOUCH)
library(PCMBase)

RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
phyltree<-ape::rtree(5)

regimes<-c("A","B","B","C","C")
regimesFitch<-fitch.mvsl(phyltree,regimes,root=1,deltran=TRUE)

## check regimes proposal
testthat::expect_equivalent(regimesFitch$branch_regimes, c("A", "A", "B", "A" ,"B", "C" ,"C", "C"))
testthat::expect_equivalent(as.character(regimesFitch$root_regime),"A")

## check if phylogeny is unchanged
testthat::expect_equivalent(regimesFitch$phyltree$edge,phyltree$edge)
testthat::expect_equivalent(regimesFitch$phyltree$tip.label,phyltree$tip.label)
testthat::expect_equivalent(regimesFitch$phyltree$Nnode,phyltree$Nnode)
testthat::expect_equivalent(regimesFitch$phyltree$edge.length,phyltree$edge.length) 
