\name{parametric.bootstrap}
\alias{parametric.bootstrap}
\title{Parametric bootstrap for confidence intervals}
\description{
The function performs a parametric bootstrap for confidence
intervals for estimates of the evolutionary model. 
The user may specify what parameters are to have 
their confidence intervals returned. The user is recommended to install 
the suggested package  \pkg{PCMBaseCpp} which significantly speeds up the 
calculations (see Details).
}
\usage{
parametric.bootstrap(estimated.model, phyltree, 
values.to.bootstrap = NULL, regimes = NULL, 
root.regime = NULL, M.error = NULL, predictors = NULL, 
kY = NULL, numboot = 100, Atype = NULL, Syytype = NULL, 
diagA = NULL, parameter_signs = NULL, start_point_for_optim = NULL,
parscale = NULL, min_bl = 0.0003, maxiter = c(10,50,100), estimateBmethod="ML")
}
\arguments{
    \item{estimated.model}{
An estimated by evolutionary model. It can be e.g. the output of \cr
\code{BrownianMotionModel()}, \code{ouchModel()}, \code{mvslouchModel()} \cr
or \code{estimate.evolutionary.model()}. 
In the last case the model under \code{BestModel} is analyzed. 
}
  \item{phyltree}{
The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a \code{nexus} file 
by the \code{read.nexus()} function from the \pkg{ape} package. The "standard" \pkg{ape} node
indexing is assumed: for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} 
and the root index \code{n+1}. The \code{root.edge} field is ignored.
}
  \item{values.to.bootstrap}{
A vector of parameter/composite statistic names that the user is interested in. They
are extracted from the bootstrapped elements for easy access.
}
  \item{regimes}{
A vector or list of regimes. If vector then each entry corresponds to each of \code{phyltree}'s
branches, i.e. to each row of \code{phyltree$edge}. If list then each list entry corresponds to a tip node
and is a vector for regimes on that lineage. If \code{NULL}, then a constant regime is assumed on
the whole tree.
}

  \item{root.regime}{
The regime at the root of the tree. If not given, then it is taken as the regime that is present
on the root's daughter lineages and is the most frequent one in the \code{regimes} vector. If more
than one regime has the same maximum frequency, then alphabetically first one of the maximum ones
is taken.
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
From version 2.0.0 of mvSLOUCH it is impossible to pass a single joint measurement error matrix 
for all the species and traits. 
}
\item{predictors}{
A vector giving the numbers of the columns from the original data
which are to be considered predictor ones, \emph{i.e.} conditioned
on in the program output. If not provided then the "X" variables 
are treated as predictors, but this only for the OUBM models (for the others in this case 
none are treated as predictors).
}
  \item{kY}{
Number of "Y" (response) variables, for the OUBM models. 
The first \code{kY}
columns of \code{mY} are the "OU" ones, while the rest the "BM" ones.
In more detail this value determines the number of columns of the (simulated) data matrix to treat as response variables ("OU" ones). For example, a 
value of 1 means that only the first column is treated as a response variable, while a value of 3 means the first three 
columns are treated as response variables. Any predictor variables ("BM" ones) the user is interested in setting for a particular model 
should therefore be placed in the final columns of the data matrix, allowing for selecting select \code{kY} columns before 
this as response variables ("OU" ones).
If \code{NULL} then
it is extracted from the provided model parameters in \code{estimated.model}.
}
  \item{numboot}{
 The number of bootstraps to perform.
}
  \item{Atype}{
The class of the \code{A} matrix. It can take one of the following values: 
\cr \code{"SingleValueDiagonal"}, \code{"Diagonal"},
\code{"UpperTri"}, \code{"LowerTri"}, \cr \code{"SymmetricPositiveDefinite"},
\code{"Symmetric"}, \code{"DecomposablePositive"}, \cr \code{"DecomposableNegative"}, 
\code{"DecomposableReal"}, \code{"Invertible"}, \code{"Any"}. If \code{NULL} then
it is extracted from the provided model parameters in \code{estimated.model}.
}
  \item{Syytype}{
The class of the Syy matrix, ignored if \code{evolmodel} equals 
\code{"BM"}. Otherwise it can take one of the following values: 
\code{"SingleValueDiagonal"}, \code{"Diagonal"},
\code{"UpperTri"}, \code{"LowerTri"},  \code{"Symmetric"}, \code{"Any"}. 
If \code{NULL} then it is extracted from the provided model parameters 
in \code{estimated.model}.
}
  \item{diagA}{
Should the diagonal of \code{A} be forced to be positive (\code{"Positive"}), \cr
negative (\code{"Negative"}) or the sign free to vary (\code{NULL}). 
However, setting this to
a non-\code{NULL} value when \code{evolmodel} is \code{"mvslouch"} might be 
(but simulations concerning this are not conclusive) slightly detrimental to the optimization process 
if \code{Atype} is \code{"DecomposablePositive"}, \code{"DecomposableNegative"}, or
\code{"DecomposableReal"}. In these cases \code{A} is parametrized by its eigendecomposition.
Additional exponentiation of the diagonal, to ensure positivity, could (but this is uncertain)
make the exploration of the likelihood surface more difficult. 
In the case of \code{Atype} being \code{"SymmetricPositiveDefinite"},
the diagonal is always guaranteed to be positive.
If \code{NULL} then the function checks if it is not in the provided model parameters 
in \code{estimated.model}.
}
   \item{parameter_signs}{
 WARNING: ONLY use this option if you understand what you are doing! This option
 is still in an experimental stage so some setups might not work (please report).
 A list allowing the user to control whether specific entries for each model parameter
 should be positive, negative, zero or set to a specific (other) value.  The entries
 of the list have to be named, the admissible names are \code{"signsA"} (for \code{A}
 matrix), \code{"signsB"} (for \code{B} matrix), \code{"signsSyy"} (for \code{Syy}
 matrix) and \code{"signsmPsi"} (for \code{mPsi} matrix) and \code{"signsvY0"} (for \code{vY0}
 matrix). Any other entry in this list will be ignored. Each entry of the list has to
 be a matrix of appropriate size, i.e. of the size of the parameter to which it corresponds.
 Inside this matrix the possible values are \code{"+"} if the given entry is to be positive,
 \code{"-"} if the given entry is to be negative, \code{x}, where \code{x} is a number,
 if the entry is to be set to specified value or \code{NA} if the entry is to be freely 
 estimated. See \code{estimate.evolutionary.model}, \code{ouchModel} and \code{mvslouchModel}
 for further details, examples and important warnings!
}
\item{start_point_for_optim}{
A named list with starting parameters for of the parameters for be optimized by \code{optim()}, currently 
only \code{A} and \code{Syy} for OUOU and OUBM models, i.e. will not work with BM model. 
One may provide both or only one of them.
Make sure that the parameter is consistent with the other parameter restrictions as no check is done and this 
can result in undefined behaviour. For example one may provide this as (provided dimensions and other
parameter restrictions agree) 
\preformatted{start_point_for_optim=list(A=rbind(c(2,0),(0,4)), 
Syy=rbind(c(1,0.5),c(0,2))).}
This starting point is always jittered in each bootstrap replicate as the employed
\code{"Nelder-Mead"} method in \code{optim()} is deterministic.
 }
  \item{parscale}{
 A vector to calculate the \code{parscale} argument for \code{optim}. It is a named
 vector with 3 entries, e.g. \cr \code{c("parscale_A"=3,"logparscale_A"=5,"logparscale_other"=1)}.
 The entry \code{parscale_A} is the scale for entries of the \code{A} matrix,
 \code{logparscale_A} is the scale for entries of the \code{A} matrix that are 
 optimized over on the logarithmic scale, e.g. if eigenvalues are assumed to be positive,
 then optimization is done over \code{log(eigenvalue)} for \code{A}'s eigendecomposition
 and \code{logparscale_other} is the scale for entries other then of \code{A}
 that are done on the logarithmic scale (e.g. \code{Syy}'s diagonal, or other entries
 indicated as positive via \code{parameter_signs}). If not provided (or if a name of the vector
 is misspelled), then made equal  to the example value provided above.  For other elements, then 
 mentioned above, that are optimized over by \code{optim()}, \code{1} is used for \code{optim()}'s 
 \code{parscale}. It is advised that  the user experiments with a couple of different values and 
 reads \code{\link[stats]{optim}}'s  man page.
 }
  \item{min_bl}{Value to which \pkg{PCMBase}'s \code{PCMBase.Threshold.Skip.Singular} should be set.
  It indicates that branches of length shorter than \code{min_bl} should be skipped in likelihood
  calculations. Short branches can result in singular covariance matrices for the transition
  density along a branch. The user should adjust this value if a lot of warnings are raised by
  \pkg{PCMBase}  about singularities during the likelihood calculations.
  Furthermore, \pkg{mvSLOUCH} sets all branches in the tree shorter than \code{min_bl} to \code{min_bl}.
  However, this does not concern tip branches-these cannot be skipped and hence should be long
  enough so that numerical issues are not raised.
  }
  \item{maxiter}{The maximum number of iterations for different components of the estimation
  algorithm. A vector of three integers. The first is the number of iterations for phylogenetic
  GLS evaluations, i.e. conditional on the other parameters, the regime optima, perhaps \code{B}, and perhaps
  initial state are estimated by a phylogenetic GLS procedure. After this the other (except of \code{B} 
  in OUBM model case) parameters are optimized over by \code{optim()}. This first
  entry controls the number of iterations of this procedure.   The second is the number of iterations inside
  the iterated GLS for the OUBM model. In the first step regime optima and \code{B} (and perhaps initial state) 
  are estimated conditional on the other parameters and current estimate of \code{B}, then the estimate of \code{B}
  is update and the same phylogenetic GLS is repeated (second entry of \code{maxiter} number of times).  
  Finally, the third is the value of \code{maxiter}  passed to \code{optim()}, apart from the optimization 
  in the Brownian motion and measurement error case. If the bootstrapped model is a Brownian motion one,
  then this parameter is ignored, if OUOU, then the second entry is ignored.
  }
  \item{estimateBmethod}{Only relevant for OUBM models, should \code{B} be estimated by maximum likelihood (default 
  value \code{"ML"}) or generalized least squares (value \code{"GLS"}).
  }
}

\details{
The likelihood calculations are done by the \pkg{PCMBase} package. However, there is a
C++ backend, \pkg{PCMBaseCpp}. If it is not available, then the likelihood is calculated
slower using pure R. However, with the calculations in C++ up to a 100-fold increase in
speed is possible (more realistically 10-20 times). The \pkg{PCMBaseCpp} package is available
from \url{https://github.com/venelin/PCMBaseCpp}.

The setting \code{Atype="Any"} means that one assumes the matrix \code{A} is eigendecomposable.
If the estimation algorithm hits a defective \code{A}, then it sets the log-likelihood at
the minimum value and will try to get out of this dip.
}

\value{
A list with all the bootstrap simulations is returned. The elements of the
list are the following.
\item{paramatric.bootstrap.estimation.replicates }{A list of length
equalling \code{numboot}. Each element is the result of the bootstrap
replicate - the estimation results in the format of the output
of \pkg{mvSLOUCH} functions, with an additional field \code{data}, 
the simulated data.}
\item{bootstrapped.parameters }{If \code{values.to.bootstrap} is not \code{NULL}
then a list of length equalling length of \code{values.to.bootstrap}. Each
element corresponds to the respective element of \code{values.to.bootstrap}
and contains a list of the bootstrapped values of this element.}
}
\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.

Butler, M.A. and A.A. King (2004) 
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.  
  American Naturalist 164:683-695.
}
\author{Krzysztof Bartoszek}

\note{
The \pkg{ouch} package implements a parametric bootstrap and reading about 
it could be helpful.
}

\section{Warning}{The estimation can take a long time and hence many bootstrap
replicates will take even more time.The code can produce (a lot of) warnings and 
errors during the search procedure, this is nothing to worry about.
}

\seealso{
\code{\link{BrownianMotionModel}}, \code{\link{estimate.evolutionary.model}},
\code{\link{mvslouchModel}}, \code{\link{ouchModel}}, \code{\link[ouch]{bootstrap}},
\code{\link[stats]{optim}}
}
\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape. 
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(5)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

BMparameters<-list(vX0=matrix(0,nrow=3,ncol=1),
Sxx=rbind(c(1,0,0),c(0.2,1,0),c(0.3,0.25,1)))

### Now simulate the data.
BMdata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
BMdata<-BMdata[phyltree$tip.label,,drop=FALSE]

### Recover the parameters of the Brownian motion.
BMestim<-BrownianMotionModel(phyltree,BMdata)

### And finally obtain bootstrap confidence intervals for some parameters
BMbootstrap<-parametric.bootstrap(estimated.model=BMestim,phyltree=phyltree,
values.to.bootstrap=c("vX0","StS"),M.error=NULL,numboot=2)
RNGversion(as.character(getRversion()))

\dontrun{##It takes too long to run this
### Define a vector of regimes.
regimes<-c("small","small","large","small","small","large","large","large")

### Define SDE parameters to be able to simulate data under the mvOUBM model.
OUBMparameters<-list(vY0=matrix(c(1,-1),ncol=1,nrow=2),A=rbind(c(9,0),c(0,5)),
B=matrix(c(2,-2),ncol=1,nrow=2),mPsi=cbind("small"=c(1,-1),"large"=c(-1,1)),
Syy=rbind(c(1,0.25),c(0,1)),vX0=matrix(0,1,1),Sxx=matrix(1,1,1),
Syx=matrix(0,ncol=1,nrow=2),Sxy=matrix(0,ncol=2,nrow=1))

### Now simulate the data.
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,regimes,NULL)
OUBMdata<-OUBMdata[phyltree$tip.label,,drop=FALSE]

### Try to recover the parameters of the mvOUBM model.
OUBMestim<-mvslouchModel(phyltree,OUBMdata,2,regimes,Atype="DecomposablePositive",
Syytype="UpperTri",diagA="Positive",maxiter=c(10,50,100))

### And finally bootstrap with particular interest in the evolutionary and optimal
### regressions

OUBMbootstrap<-parametric.bootstrap(estimated.model=OUBMestim,phyltree=phyltree,
values.to.bootstrap=c("evolutionary.regression","optimal.regression"),
regimes=regimes,root.regime="small",M.error=NULL,predictors=c(3),kY=2,
numboot=5,Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive",
maxiter=c(10,50,100))


### We now demonstrate an alternative setup
### Define SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),"large"=c(-1,1,0.5)),
Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

### Now simulate the data.
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

### Try to recover the parameters of the OUOU model.
estimResults<-estimate.evolutionary.model(phyltree,OUOUdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=3,model.setups=NULL,predictors=c(3),kY=2,
doPrint=TRUE,pESS=NULL,maxiter=c(10,50,100))

### And finally bootstrap with particular interest in the evolutionary regression
OUOUbootstrap<-parametric.bootstrap(estimated.model=estimResults,phyltree=phyltree,
values.to.bootstrap=c("evolutionary.regression"),
regimes=regimes,root.regime="small",M.error=NULL,predictors=c(3),kY=NULL,
numboot=5,Atype=NULL,Syytype=NULL,diagA=NULL)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{multivariate}
\keyword{htest}
