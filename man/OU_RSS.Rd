\name{OU_RSS}
\alias{OU_RSS}
\title{
Calculates the RSS under a (multivariate) phylogenetic OU model of evolution
}
\description{
The \code{OU_RSS} function calculates the residual sum of squares (RSS)
for given data under a multivariate OU model evolving on the phylogeny. 
The user is recommended to install the suggested package 
\pkg{PCMBaseCpp} which significantly speeds up the calculations (see Details).
}
\usage{
OU_RSS(mY, phyltree, modelParams, M.error = NULL, do_centre = NA,
regimes = NULL, regimes.times = NULL, root.regime = NULL)
}
\arguments{
  \item{mY}{
A matrix with the rows corresponding to the tip species while the columns correspond to the traits.
The rows should be named by species \cr (field \code{phyltree$tip.label}), if not, then a warning
is thrown and the order of the species is assumed to be the same as the order in which the 
species are in the phylogeny (i.e. correspond to the node indices \code{1:n}, 
where \eqn{n}{n} is the number of tips). The columns should be named by traits, otherwise
a warning is thrown and generic names are generated.
}
  \item{phyltree}{
The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a \code{nexus} file 
by the \code{read.nexus()} function from the \pkg{ape} package. The "standard" \pkg{ape} node
indexing is assumed: for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} 
and the root index \code{n+1}. The \code{root.edge} field is ignored.
}
  \item{modelParams}{
List of model parameters of the BM/OUOU/OUBM model as \code{ParamsInModel} part of output of 
\code{BrownianMotionModel}/\code{ouchModel}/\code{mvslouchModel}.
}
  \item{M.error}{
An optional measurement error covariance structure. The measurement errors between species are assumed 
independent. The program tries to recognize the structure of the passed matrix and accepts the following 
possibilities :
\itemize{
\item{}{a single number that is a common measurement error for all tips and species,}
\item{}{a m element vector with each value corresponding to a variable, measurement errors are independent
between variables and each species is assumed to have the same measurement errors,}
\item{}{a m x m ((number of variables) x (number of variables)) matrix, all species will
have the same measurement error,}
\item{}{a list of length n (number of species), each list element is the covariance structure 
for the appropriate (numbering according to tree) species, either a single number 
(each variable has same variance), vector (of length m for each variable), or m x m matrix, 
the order of the list has to correspond to the order of the nodes in the \code{phyltree} object,}
\item{}{NULL no measurement error.}
}
From version \code{2.0.0} of \pkg{mvSLOUCH} it is impossible to pass a single joint measurement error matrix 
for all the species and traits.
}
  \item{do_centre}{
Should the data, \code{mY}, be centred (admissable values are \code{"average"}, \code{"phylaverage"} or 
\code{"evolutionary_model"}) or not (\code{NA}). If \code{"average"}, then each column (trait)
is centred by its arithmetic average, if \code{"evolutionary_model"}, then \code{mY} is centred
by the expectation under the evolutoniary model (inferred from \code{modelParams}) and if
\code{"phylaverage"}, then \code{mY} is centred by a phylogenetically weighted arithmetic
average (see Details).
}

  \item{regimes}{
A vector or list of regimes. If vector then each entry corresponds to each of the branches of \code{phyltree},
i.e. to each row of \code{phyltree$edge}. If list then each list entry corresponds to a tip node
and is a vector for regimes on that lineage. If \code{NULL}, then a constant regime is assumed on
the whole tree.
}
  \item{regimes.times}{
A list of vectors for each tree node, it starts with \eqn{0}{0} and ends with the current time of the species.
In between are the times where the regimes (niches) changed. If \code{NULL} then each branch is
considered to be a regime.
}
  \item{root.regime}{
The regime at the root of the tree. If not given, then it is taken as the regime that is present
on the daughter lineages stemming from the root and is the most frequent one in the \code{regimes} 
vector. If more than one regime has the same maximum frequency, then alphabetically first one of 
the maximum ones is taken.
}
}

\details{
The matrix algebra calculations are done using the likelihood function offerred by the \pkg{PCMBase} package. 
However, there is a C++ backend, \pkg{PCMBaseCpp}. If it is not available, then the likelihood is calculated
slower using pure R. However, with the calculations in C++ up to a 100-fold increase in
speed is possible (more realistically 10-20 times). The \pkg{PCMBaseCpp} package is available
from \url{https://github.com/venelin/PCMBaseCpp}.

For a given input data matrix, \code{mY}, the function considers the stacking of it by rows
(i.e. stacking species by species). Let Y = vec(\code{mY}), i.e. \code{Y<-c(t(mY))},
V be the between-species-between-traits variance-covariance matrix (under the parameters passed
in \code{modelParams}) and \code{v} a centring vector (if \code{do_centre} is \code{NA}, then
\eqn{v=0}{v=0}). The function calculates the value of the quadratic form
\deqn{(Y-v)^{T}V^{-1}(Y-v).}{\code{t(Y-v)\%*\%solve(V)\%*\%(Y-v)}.}

A special centring is when \code{do_centre} equals \code{"phylaverage"}. In this situation
the centring vector is a phylogenetically weighted average, i.e.
\deqn{v=(D^{T}V^{-1}D)^{-1}D^{T}V^{-1}Y ,}{\code{v<-solve(t(D)\%*\%solve(V)\%*\%D)\%*\%t(D)\%*\%solve(V)\%*\%Y}, }
where denoting \eqn{1_{n}}{1n} as a column vector of n ones and \eqn{Id_{k}}{diag(k)} as 
the identity matrix with rows and columns equalling the number or columns of \code{mY},
\deqn{D=1_{n}\otimes Id_{k}.}{\code{D<-rep(1,n)\%x\%diag(k)}.}
}
\value{
The value of the residual sum of squares, quadratic form with respect to the 
between-species-between-traits precision matrix. Also the used phylogeny is returned.
}
\author{Krzysztof Bartoszek}

\examples{
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
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

## below will return the RSS under the assumed OUOU model of evolution
OU_RSS(OUOUdata, phyltree, OUOUparameters, M.error=NULL, 
do_centre="evolutionary_model", regimes = regimes)
}

\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{array}
\keyword{algebra}
