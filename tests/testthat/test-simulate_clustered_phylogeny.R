## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(testthat)
context("mvSLOUCH: simulate_clustered_phylogeny")

library(mvSLOUCH)
library(PCMBase)

RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
## We use a wrapper function for illustration
## a single phylo object
my_sim.bd.taxa<-function(n,...){
    ape::rphylo(n=n,...)
}
phyltree1<-simulate_clustered_phylogeny(v_sizeclusts=c(5,5,5),f_simclustphyl=my_sim.bd.taxa,
b_change_joining_branches=TRUE, joining_branchlengths=c(20,NA),joining=my_sim.bd.taxa,
birth=1,death=0)

RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
## The below code should return the same tree as above
phyltree2<-simulate_clustered_phylogeny(v_sizeclusts=c(5,5,5),f_simclustphyl="sim.bd.taxa_Yule1",
b_change_joining_branches=TRUE, joining_branchlengths=c(20,NA),joining="sim.bd.taxa_Yule1")

testthat::expect_equivalent(phyltree1,phyltree2)
testthat::expect_identical(length(phyltree1$tips_clusters),3L)
testthat::expect_identical(length(phyltree1$edges_clusters),length(phyltree1$tips_clusters)+1L)
testthat::expect_identical(as.integer(phyltree1$Nnode)+1L,length(phyltree1$tip.label))
testthat::expect_identical(length(phyltree1$tip.label),15L)

