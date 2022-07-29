\name{OU_phylreg}
\alias{OU_phylreg}
\title{
Performs a phylogenetic regression under a given OU model of evolution
}
\description{
The \code{OU_phylreg} function does a phylogenetic regression
for given response and design matrices under a multivariate OU model evolving on the phylogeny. 
The user is recommended to install the suggested package 
\pkg{PCMBaseCpp} which significantly speeds up the calculations (see Details).
}
\usage{
OU_phylreg(mY, mD, phyltree, modelParams, regimes = NULL, kY = NULL, M.error = NULL,
signif_level = 0.05, regimes.times = NULL, root.regime = NULL, b_GLSB = FALSE, 
b_GLSX0 = FALSE, signsB = NULL, signsvX0 = NULL, estimate.root.state = FALSE)
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
  \item{mD}{
A design matrix with the rows corresponding to the traits in the tips species while the columns correspond to 
the unknown regression variables. The number or rows have to correspond to the number of elements in
\code{mY}, the data are assumed to be stacked by species. If \code{NA} it is assumed to be the 
design matrix to estimate regression parameters under the given model of evolution, see Details.
If it is \code{"phylaverage"}, then a phylogenetically weighted average is calculated, see Details.
}

  \item{phyltree}{
The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a \code{nexus} file 
by the \code{read.nexus()} function from the \pkg{ape} package. The "standard" \pkg{ape} node
indexing is assumed: for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} 
and the root index \code{n+1}. The \code{root.edge} field is ignored.
}
  \item{modelParams}{
List of model parameters of the BM/OUOU/OUBM model as \code{ParamsInModel} part of output of 
\code{BrownianMotionModel}/\code{ouchModel}/\code{mvslouchModel}. Some of them can be \code{NA}
in order to be estimated by the regression procedure, see Details.
}
  \item{regimes}{
A vector or list of regimes. If vector then each entry corresponds to each of the branches of \code{phyltree},
i.e. to each row of \code{phyltree$edge}. If list then each list entry corresponds to a tip node
and is a vector for regimes on that lineage. If \code{NULL}, then a constant regime is assumed on
the whole tree.
}

  \item{kY}{Number of "Y" (response) variables if the considered model is an OUBM one. The first \code{kY}
columns of \code{mY} are the "OU" ones, while the rest the "BM" ones.
In more detail this value determines the number of columns of the \code{mY} matrix to treat as response variables ("OU" ones). For example, a 
value of 1 means that only the first column is treated as a response variable, while a value of 3 means the first three 
columns are treated as response variables. Any predictor variables ("BM" ones) the user is interested in setting for a particular model 
should therefore be placed in the final columns of the \code{mY} matrix, allowing for selecting select \code{kY} columns before 
this as response variables ("OU" ones).
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

  \item{signif_level}{
The significance level to be taken when calculating regression confidence intervals, i.e. 
\eqn{(1-signif_level)\cdot 100}{(1-signif_level)*100} percent confidence intervals are returned.
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

  \item{b_GLSB}{
If the evolutionary model is an OUBM one (and \code{mD} is \code{NA}), then should the \eqn{B}{\code{B}}
matrix be also estimated. If \eqn{B}{\code{B}} is completely provided in \code{modelParams}, then that 
value is taken as an initial guess for the regression estimation procedure (as the 
between-species-between-traits variance-covariance matrix depends on \code{B}).
}

  \item{b_GLSX0}{
If the evolutionary model is an OUBM or BM one (and \code{mD} is \code{NA}), then should the \eqn{X_{0}}{\code{vX0}}
ancestral vector  also be estimated. If \eqn{X_{0}}{\code{vX)}} is completely provided in \code{modelParams}, then that 
value is taken as an initial guess for the regression estimation procedure (as in the OUBM the design matrix
depends on \code{vX0}, in the BM case this value is ignored).
}

  \item{signsB}{
A matrix of constraints on the estimation of \code{B}, with the same dimensions as \code{B}, 
if \code{b_GLBS} is \code{TRUE}.  Inside this matrix the possible values are \code{"+"} 
if the given entry is to be positive,  \code{"-"} if the given entry is to be negative, 
\code{x}, where \code{x} is a number, if the entry is to be set to specified value or \code{NA} if the entry is 
to be freely estimated.  This option  is still in an experimental stage so some setups might not work (please report).
}

  \item{signsvX0}{
A matrix of constraints on the estimation of \code{vX0}, with the same dimensions as \code{vX0}, 
if \code{b_GLBX0} is \code{TRUE}.  Inside this matrix the possible values are \code{"+"} 
if the given entry is to be positive,  \code{"-"} if the given entry is to be negative, 
\code{x}, where \code{x} is a number, if the entry is to be set to specified value or \code{NA} if the entry is 
to be freely estimated.  This option  is still in an experimental stage so some setups might not work (please report).
}

\item{estimate.root.state}{
    Should the root state be estimate \code{TRUE} (not recommended) or set at the optimum 
    \code{FALSE} (recommended). Root state estimation is usually unreliable hence
    if fossil measurements are available prediction based on them and the estimated
    model will probably be more accurate. If there is only one regime, then 
    estimation of the root state separately is impossible and will not be allowed.
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
in \code{modelParams}). The function calculates the value of the generalized least squares 
estimator (not directly, but as a transformation of the likelihood provided by \pkg{PCMBase})
\deqn{v=(D^{T}V^{-1}D)^{-1}D^{T}V^{-1}Y.}{\code{v<-solve(t(D)\%*\%solve(V)\%*\%D)\%*\%t(D)\%*\%solve(V)\%*\%Y}.}

The user can provide the design matrix directly or if \code{mD} is \code{NA}, then the design matrix
induced by the evolutionary model in \code{modelParams} is assumed. The following parameters
can be estimated: \code{vX0} (if \code{b_GLSX0} is \code{TRUE}, BM model); \code{mPsi}, 
\code{vY0} (if \code{estimate.root.state} is \code{TRUE}, otherwise set at optimum) for OUOU model; 
\code{vX0} (if \code{b_GLSX0} is \code{TRUE}) \code{mPsi}, \code{vY0} (if \code{estimate.root.state} 
is \code{TRUE}, otherwise set at optimum), \code{B} (if \code{b_GLSB} is \code{TRUE}). One can
constrain (some of) the elements of the matrices to be estimated to be postive, negative or 
equal to some value. For \code{B} and \code{vX0} this was described in the description of the 
arguments of \code{signsB} and \code{signsvX0}. For \code{mPsi} and \code{vY0} one does this
in the respective entries of \code{modelParams}. There matrix entries can be set to \code{"+"},
\code{"-"}, \code{NA} or some specific value. In the OUBM case the model specfic design matrix
is not derived from the conditional expectation of all of the responses on all of the predictors,
but from the conditional expectations of each tip species independently (as if V were block
diagonal). This is as the joint condtional expecation design matrix cannot be calculated 
at the moment in an efficient manner and would cause a serious computational bottleneck.
However this only makes a difference if \code{B} is to be estimated inside the GLS. 

Special support is given if one wants to compute a phylogenetically weighted mean.
If \code{mD} is set to \code{"phylaverage"}, then it is calculated as 
\deqn{D_p=1_{n}\otimes Id_{k}, }{\code{Dp<-rep(1,n)\%x\%diag(k)}, }
where \eqn{1_{n}}{1n} is a column vector of n ones and \eqn{Id_{k}}{diag(k)} is
the identity matrix with rows and columns equalling the number or columns of \code{mY}.
}

\value{
A list with the following entries
\itemize{
\item{vGLSest }{The regression estimates}
\item{regression.covariance.matrix }{The covariance matrix between regression estimates.}
\item{regression.confidence.intervals }{The confidence intervals for each estimated parameter.}
\item{modelParams }{The model parameters updated if anything was estimated from them in the procedure.}
\item{mD }{The used or calculated design matrix.}
\item{RSS }{The residual sum of squares.}
\item{R2_average }{R2, where the alternative model is the sample average.}
\item{R2_phylaverage }{R2, where the alternative model is the phylogenetically weighted sample average, i.e. the design matrix
is \eqn{D_p }{\code{Dp}}.}
\item{RSS_average }{The RSS with respect to the sample average.}
\item{RSS_phylaverage }{The RSS with respect to the phylogenetically weighted sample average.}
\item{phyltree }{The phylogeny used, returned as in the estimation procedure some additional fields are calculated. This could 
help in a speed up if the \code{OU_phylreg} is used in some iterative procedure.}
}
}

\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.
  
Hansen, T.F. and Bartoszek, K. (2012)
  Interpreting the evolutionary regression: the interplay between observational and 
  biological errors in phylogenetic comparative studies.
  Systematic Biology 61(3):413-425.
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

OUOUparameters_reg<-OUOUparameters
OUOUparameters_reg$mPsi<-apply(OUOUparameters_reg$mPsi,c(1,2),function(x){NA})
OUOUparameters_reg$vY0<-apply(OUOUparameters_reg$vY0,c(1,2),function(x){NA})
## estimate parameters under OUOU model
OU_phylreg(OUOUdata, NA, phyltree, OUOUparameters_reg, regimes=regimes, 
kY=NULL, M.error=NULL)

}
\keyword{models}
\keyword{regression}
\keyword{multivariate}
\keyword{array}
\keyword{algebra}
