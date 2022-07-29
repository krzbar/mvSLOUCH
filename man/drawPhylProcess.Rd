\name{drawPhylProcess}
\alias{drawPhylProcess}

\title{Plots the realization of a process evolving on a phylogenetic tree}

\description{
The function takes the output of the simulation functions
and based on it plots the realization of the process on the tree.
Can handle multiple traits, in this case each trait is plotted separately.
The function does draw anything else (like axes) but the realization
of the process. Any additions are up to the user.
}
\usage{
drawPhylProcess(PhylTraitProcess, vTraitsToPlot=NULL, vColours = "black", 
plotlayout = c(1, 1), additionalfigs = FALSE, modelParams = NULL, 
EvolModel = NULL, xlimits = NULL, ylimits = NULL)
}

\arguments{
  \item{PhylTraitProcess}{
The simulated realization of the process, the direct output of one of the package's simulation 
function or a matrix (if \code{fullTrajectory} is \code{TRUE}).
In the second case the matrix consists of k+1 columns, where k is the number of traits. The first column
are the time instances, the next k the values of the traits at that instance. Since evolution
takes place on a phylogenetic tree - there should be multiple copies of the same time moment, i.e.
one for each branch of the tree.
}
  \item{vTraitsToPlot}{A vector providing the column numbers of the traits to plot. If \code{NULL}, then all
  traits are plotted. The column numbers have to be obtained from the \cr \code{PhylTraitProcess} object,
  the matrix \code{$trajectory} for each branch. Notice that the first column is time! The same trait
  may be plotted multiple times (but a warning will be raised).}
  \item{vColours}{A vector of colours to be used for each trait. If length is less than the number of 
  traits then colours are recycled}
  \item{plotlayout}{How many plots per page if more than one trait, i.e. \code{par(mfrow=plotlayout))}.}
  \item{additionalfigs}{Should additional items be plotted on each figure, the ancestral state and 
  deterministic, 
  \eqn{\Psi} when appropriate. If there are many regime levels then only the first column of 
  \eqn{\Psi} is used.}
  \item{modelParams}{List of model parameters.}
  \item{EvolModel}{The evolutionary model.}
  \item{xlimits}{The x limits of the plot. Can be useful to fix if one wants to have a number of graphs
  on the same scale. This can be either a vector of length 2 (minimum and maximum value of the x-axis),
  or a list of length equalling the number of traits with each entry being a vector of length 2 or 
  a matrix with two columns and rows equalling the number of traits. If not provided then the value
  is just the minimum and maximum from the data.} 
  \item{ylimits}{The y limits of the plot. Can be useful to fix if one wants to have a number of graphs
  on the same scale. This can be either a vector of length 2 (minimum and maximum value of the x-axis),
  or a list of length equalling the number of traits with each entry being a vector of length 2 or 
  a matrix with two columns and rows equalling the number of traits. If not provided then the value
  is just the minimum and maximum from the data.} 
}

\value{Returns a meaningless NA value.}
\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.
}

\author{Krzysztof Bartoszek}


\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape.
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(3)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

### Define a vector of regimes.
#regimes<-c("small","small","large","small","small","large","large","large")
#regimes<-c("small","small","large","small","small","large")
regimes<-c("small","small","large","small")

### Define SDE parameters to be able to simulate data under the OUOU model.
## 3D model
## OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
## A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),
## "large"=c(-1,1,0.5)),Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

## 2D model for speed on CRAN
OUOUparameters<-list(vY0=matrix(c(1,-1),nrow=2,ncol=1),
A=rbind(c(9,0),c(0,5)),mPsi=cbind("small"=c(1,-1),
"large"=c(-1,1)),Syy=rbind(c(1,0.25),c(0,1)))


### Now simulate the data keeping the whole trajectory
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL,fullTrajectory=TRUE)
drawPhylProcess(PhylTraitProcess=OUOUdata,plotlayout=c(1,3))
RNGversion(as.character(getRversion()))

}

\keyword{hplot}
