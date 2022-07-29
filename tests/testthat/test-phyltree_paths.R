## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: phyltree_paths")

library(mvSLOUCH)
library(PCMBase)

RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
phyltree<-ape::rtree(5)
phyltree_augmented<-phyltree_paths(phyltree)

testthat::expect_identical(length(phyltree_augmented$path.from.root),phyltree_augmented$Ntips)
testthat::expect_identical(length(phyltree_augmented$time.of.nodes),phyltree_augmented$Nnode+phyltree_augmented$Ntips)

## check if phylogeny is unchanged
testthat::expect_equivalent(phyltree_augmented$edge,phyltree$edge)
testthat::expect_equivalent(phyltree_augmented$tip.label,phyltree$tip.label)
testthat::expect_equivalent(phyltree_augmented$Nnode,phyltree$Nnode)
testthat::expect_equivalent(phyltree_augmented$edge.length,phyltree$edge.length) 
