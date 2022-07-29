\name{plot.clustered_phylo}
\alias{plot.clustered_phylo}

\title{Plots a clustered_phylo object.}

\description{
The function plots a clustered_phylo object allowing the user to differently
visualize the different clades/clusters on the phylogeny and also the
joining them subtree.
}
\usage{
\method{plot}{clustered_phylo}(x, clust_cols = NULL, clust_edge.width = NULL, 
clust_edge.lty = NULL, clust_tip.color = "black", joiningphylo_col = "black", 
joiningphylo_edge.width = 1, joiningphylo_edge.lty = 1, ...)
}

\arguments{
    \item{x}{A phylogenetic tree of class \code{clustered_phylo}, i.e. output of 
    \code{mvSLOUCH::simulate_clustered_phylogeny()}}
    \item{clust_cols}{Vector of colours of edges inside each cluster. Default \code{NULL}, corresponding to
    \code{"black"}. If length of this vector does not equal to the number of clusters, then it is recycled.}
    \item{clust_edge.width}{Numeric vector of widths of edges inside each cluster. Default \code{NULL}, corresponding to
    \code{1}. If length of this vector does not equal to the number of clusters, then it is recycled.}
    \item{clust_edge.lty}{Vector of an edge's type inside each cluster. Default \code{NULL}, corresponding to
    \code{1}. If length of this vector does not equal to the number of clusters, then it is recycled.}
    \item{clust_tip.color}{Vector of colours of tips' labels inside each cluster. Default \code{NULL}, corresponding to
    \code{"black"}. If length of this vector does not equal to the number of clusters, then it is recycled.}
    \item{joiningphylo_col}{Colour of edges inside the subtree joining the clusters.} 
    \item{joiningphylo_edge.width}{Width of edges inside the subtree joining the clusters.}
    \item{joiningphylo_edge.lty}{Edges' type inside the subtree joining the clusters.}
    \item{...}{Other parameters to be passed to \code{plot.phylo()}. Notice that here we cannot have
    \code{edge.color}, \code{edge.width} and \code{edge.lty}.}
}

\value{Same as \code{plot.phylo()}.}

\author{Krzysztof Bartoszek}


\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")

phyltree<-simulate_clustered_phylogeny(v_sizeclusts=c(5,5,5),f_simclustphyl="sim.bd.taxa_Yule1",
b_change_joining_branches=TRUE, joining_branchlengths=c(20,NA),joining="sim.bd.taxa_Yule1")

plot(phyltree,clust_cols=c("red","green","blue"),clust_edge.width=3,clust_edge.lty=c(1,2,3),
clust_tip.color=c("red","blue","green"),joiningphylo_col="black",joiningphylo_edge.width=3,
joiningphylo_edge.lty=1)

## and not plot without tip labels
plot(phyltree,clust_cols=c("red","green","blue"),clust_edge.width=3,clust_edge.lty=c(1,2,3),
joiningphylo_col="black",joiningphylo_edge.width=3,joiningphylo_edge.lty=1,show.tip.label=FALSE)

RNGversion(as.character(getRversion()))
}

\keyword{hplot}
