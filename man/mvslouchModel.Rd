\name{mvslouchModel}
\alias{mvslouchModel}
\title{
Estimate parameters under a (multivariate) OUBM model of evolution
}
\description{
The \code{mvslouchModel} function uses maximum likelihood to fit parameters of a multivariate 
OUBM model evolving on the phylogeny. The user is recommended to install the suggested package 
\pkg{PCMBaseCpp} which significantly speeds up the calculations (see Details).
}
\usage{
mvslouchModel(phyltree, mData, kY, regimes = NULL, regimes.times = NULL, 
root.regime = NULL, predictors = NULL, M.error = NULL, Atype = "Invertible", 
Syytype = "UpperTri", diagA = "Positive", estimate.root.state=FALSE, 
parameter_signs=NULL, start_point_for_optim = NULL, parscale = NULL, 
min_bl = 0.0003, maxiter = c(10,50,100), estimateBmethod="ML")
}
\arguments{
  \item{phyltree}{
The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a \code{nexus} file 
by the \code{read.nexus()} function from the \pkg{ape} package. The "standard" \pkg{ape} node
indexing is assumed: for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} 
and the root index \code{n+1}. The \code{root.edge} field is ignored.
}
  \item{mData}{
A matrix with the rows corresponding to the tip species while the columns correspond to the traits.
The rows should be named by species \cr (field \code{phyltree$tip.label}), if not, then a warning
is thrown and the order of the species is assumed to be the same as the order in which the 
species are in the phylogeny (i.e. correspond to the node indices \code{1:n}, 
where \eqn{n}{n} is the number of tips). The columns should be named by traits, otherwise
a warning is thrown and generic names are generated.
}
  \item{kY}{Number of "Y" (response) variables.
The first \code{kY}
columns of \code{mY} are the "OU" ones, while the rest the "BM" ones.
In more detail this value determines the number of columns of the \code{mData} matrix to treat as response variables ("OU" ones). For example, a 
value of 1 means that only the first column is treated as a response variable, while a value of 3 means the first three 
columns are treated as response variables. Any predictor variables ("BM" ones) the user is interested in setting for a particular model 
should therefore be placed in the final columns of the \code{mData} matrix, allowing for selecting select \code{kY} columns before 
this as response variables ("OU" ones).
 }
  \item{regimes}{
A vector or list of regimes. If vector then each entry corresponds to each of \code{phyltree}'s
branches, i.e. to each row of \code{phyltree$edge}. If list then each list entry corresponds to a tip node
and is a vector for regimes on that lineage. If \code{NULL}, then a constant regime is assumed on
the whole tree.
}
  \item{regimes.times}{
A list of vectors for each tree node, it starts with 0 and ends with the current time 
of the species. In between are the times where the regimes (niches) changed. If \code{NULL} 
then each branch is considered to be a regime.
}
  \item{root.regime}{
The regime at the root of the tree. If not given, then it is taken as the regime that is present
on the root's daughter lineages and is the most frequent one in the \code{regimes} vector. If more
than one regime has the same maximum frequency, then alphabetically first one of the maximum ones
is taken.
}
  \item{predictors}{
A vector giving the numbers of the columns from \code{mData} matrix which are to be considered predictor ones, 
\emph{i.e.} conditioned on in the program output. If not provided then columns \code{(kY+1):ncol(mData)}, i.e. the 
"BM" ones, are treated as predictors.
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
  \item{Atype}{
What class does the A matrix in the multivariate OUBM model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"},
\code{"Symmetric"}, \code{"SymmetricPositiveDefinite"}, \cr
\code{"DecomposablePositive"},  \code{"DecomposableNegative"}, \cr
\code{"DecomposableReal"}, \code{"Invertible"}, \code{"TwoByTwo"}, \code{"Any"}
}
  \item{Syytype}{
What class does the Syy matrix in the multivariate OUBM model belong to, possible values :
\code{"SingleValueDiagonal"}, \code{"Diagonal"}, \code{"UpperTri"}, \code{"LowerTri"},
\code{"Symmetric"}, \code{"Any"}
}
   \item{diagA}{
  Whether the values on \code{A}'s diagonal are to be \code{"Positive"}, \code{"Negative"} 
  or sign allowed to vary, \code{NULL}. However, setting this to
a non-\code{NULL} value might be (but simulations concerning this are not conclusive) 
slightly detrimental to the optimization process if \code{Atype} is
\code{"DecomposablePositive"}, \code{"DecomposableNegative"}, or
\code{"DecomposableReal"}. In these cases \code{A} is parametrized by its eigendecomposition.
Additional exponentiation of the diagonal, to ensure positivity, could (but this is uncertain)
make the exploration of the likelihood surface more difficult. The user is advised to
also try \code{diag=NULL}. In the case of \code{Atype} being \code{"SymmetricPositiveDefinite"},
the diagonal is always guaranteed to be positive.
}
\item{estimate.root.state}{
    Should the root state be estimate \code{TRUE} (not recommended) or set at the optimum 
    \code{FALSE} (recommended). Root state estimation is usually unreliable hence
    if fossil measurements are available prediction based on them and the estimated
    model will probably be more accurate. If there is only one regime, then 
    estimation of the root state separately is impossible and will not be allowed.
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
 estimated. See Details for an example, further description and important warnings! 
}
 \item{start_point_for_optim}{
  A name list with starting parameters for of the parameters for be optimized by \code{optim()},
  in this case \code{A} and \code{Syy}. One may provide both or only one of them.
  Make sure that the parameter is consistent with the other parameter restrictions
  as no check is done and this can result in undefined behaviour.
 For example one may provide this as (provided dimensions and other parameter restrictions agree) \preformatted{
start_point_for_optim=list(A=rbind(c(2,0),(0,4)), 
Syy=rbind(c(1,0.5),c(0,2))).}
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
  However, this does not concern tip branches-these cannot be skipped and hence should be long
  enough so that numerical issues are not raised.
  }
   \item{maxiter}{The maximum number of iterations for different components of the estimation
  algorithm. A vector of three integers. The first is the number of iterations for phylogenetic
  GLS evaluations, i.e. conditional on the other parameters, the regime optima, \code{B} and perhaps
  initial state are estimated by a phylogenetic GLS procedure. After this the other (except of \code{B}) 
  parameters are optimized over by \code{optim()}. This first entry controls the number of iterations of this 
  procedure.  The second is the number of iterations inside the iterated GLS. In the first step regime optima 
  and \code{B} (and perhaps initial state) are estimated conditional on the other parameters and current estimate 
  of \code{B}, then the estimate of \code{B} is update and the same phylogenetic GLS is repeated 
  (second entry of \code{maxiter} number of times). Finally, the third is the value of \code{maxiter}  
  passed to \code{optim()}, apart from the optimization in the Brownian motion and measurement error case.
  }
  \item{estimateBmethod}{Should \code{B} be estimated by maximum likelihood (default 
  value \code{"ML"}) or generalized least squares (value \code{"GLS"}).
  }
}
\details{
The likelihood calculations are done by the \pkg{PCMBase} package. However, there is a
C++ backend, \pkg{PCMBaseCpp}. If it is not available, then the likelihood is calculated
slower using pure R. However, with the calculations in C++ up to a 100-fold increase in
speed is possible (more realistically 10-20 times). The \pkg{PCMBaseCpp} package is available
from \url{https://github.com/venelin/PCMBaseCpp}.

This function estimates the parameters of the following multivariate SDE,
\deqn{\begin{array}{rclccl}
dY(t) & = & -A(Y(t)-(\Psi(t)- A^{-1}BX(t)))dt + \Sigma_{yy} dB(t) & Y(0) & = & Y_{0}, 
\\ dX(t) & = & \Sigma_{xx} dB(t) & X(0) & = & X_{0}
\end{array}
}{
dY(t) = -A(Y(t)-(\Psi(t)- solve(A)BX(t)))dt + Syy dB(t) Y(0)=Y0, 

dX(t) = Sxx dB(t), X(0)=X0
}
on a phylogenetic tree. It uses a numerical optimization over \code{A} (parametrized by its 
eigenvalues and eigenvectors or its QR decomposition) and \code{S} (parametrized by its values) 
and conditional on \code{A} and \code{S} estimates the values of \code{Psi} corresponding
to the different regimes by a GLS estimate. \code{Y(0)} is assumed to be equal to 
\code{- solve(A)BX(0)} plus the root value of \code{Psi}. This assumes that \code{A} is invertible. 
If not, then \code{Y(0)} will be set at the root value of \code{Psi}. This is unless 
\code{estimate.root.state=TRUE}, in such a case \code{Y(0)} will be estimated by least squares.

The setting \code{Atype="Any"} means that one assumes the matrix \code{A} is eigendecomposable.
If the estimation algorithm hits a defective \code{A}, then it sets the log-likelihood at
the minimum value and will try to get out of this dip.

The function parameter \code{parameter_signs} is special in the sense that it can give 
the user great control over the estimation procedure but can also make the output
very inconsistent with what the user provides. If we have two response traits (OU ones)
and two predictor traits (BM ones), then an EXAMPLE setting of this can be:
\cr
\code{parameter_signs=list(signsA=rbind(c("+","-"),c(0,"+")),} \cr
\code{signsSyy=rbind(c(NA,0),c(0,NA)), signsB=rbind(c(NA,0),c(0,NA)))}. 
This means that \code{A} is upper triangular with positive
values on the diagonal and a negative value on the off-diagonal, \code{Syy} is diagonal
and \code{B} is also diagonal. It is advisable to set now \code{Atype="Any"} and
\code{Syytype="Any"} (see further description).

If the given model parameter is to be estimated
by a generalized least squares (currently \code{B}, \code{mPsi} and \code{vY0}), then the 
sign specifications are ignored. However, it is possible to set specific values.
Furthermore, the package does not check (for \code{A}
and \code{Syy}) if the specifications here agree with the \code{Atype}, \code{Syytype}
and \code{diagA}. The settings in \code{signsA} and \code{signsSyy} will override
the other settings. Hence, it is up to the user to make sure that the settings of 
\code{signsA} and \code{signsSyy} are consistent with \code{Atype}, \code{Syytype}
and \code{diagA}. It is advisable to use \code{signsA} with \code{"+"}
on the diagonal and have \code{diagA=NULL}. The diagonal of \code{Syy} is forced to
be positive (unless \code{"-"} is used on the diagonal of \code{signsSyy}
but this is strongly discouraged) so it is advisable to keep \code{NA} 
on the diagonal of \code{signsSyy} and not put there \code{"+"} there.
Hence, in particular using the signs mechanism result in a wrong class of the matrix
\cr 
(e.g. \code{Atype="SymmetricPositiveDefinite"}, but after corrections for the provided entries in
\code{signsA} one obtains a non-symmetric \code{A} with complex, negative-real-part eigenvalues). 
Lastly, using \code{signsA} and \code{signsSyy} can result in
a wrong amount of \code{dof} and in turn incorrect \code{AICc} and \code{BIC} values.
What the code does is subtracts the amount of fixed values in \code{signsA} and \code{signsSyy}
from the amount of free parameters used to estimate \code{A} and \code{Syy}. For example
if one sets \cr \code{Atype="SingleValueDiagonal"} (estimated by one free parameter) 
but specified two off-diagonal values, then the amount of dofs from \code{A} will be \code{-1}!!
The ONLY fail-safe way to use this is to set \code{Atype="Any"} (if \code{signsA} used) and 
\code{Syytype="Any"} (if \code{signsSyy} used). If using \code{Syytype="Any"} and \code{signsSyy} 
the it is strongly advisable to set the entries either below or above \code{Syy}'s diagonal to \code{0}.
The reason is that \eqn{\Sigma_{yy}\Sigma_{yy}^{T}}{Syy t(Syy)} enters the likelihood and not the 
given value of \eqn{\Sigma_{yy}}{Syy}. Hence, having values below (or respectively above) the diagonal
results in an overparameterized model. The package has the option of mixing different matrix types
with specifying values in it but this is only for advanced users who need to dig into the code
to see what the \code{dof}'s should be and  if it is possible to find a correspondence between the 
parametrization and settings. If entries of \code{mPsi}, \code{vY0} and \code{B} are pre-specified,
then the \code{dof} are correctly adjusted for this. The estimation procedures currently ignore any
pre-specified values for \code{vX0} and \code{Sxx}!


The found point is described by a list containing four fields. 
The first field \cr \code{HeuristicSearchPointFinalFind} is the parametrization of the model 
parameters at the considered
point with the value of the log--likelihood.
The field \code{ParamsInModel} is the point estimate of the parameters of the SDE. 
The field \code{ParamSummary} are different composite (evaluated at the tree's height) and summary statistics,
The field \code{phylhalflife} are the eigenvalues, eigenvectors and phylogenetic half lives
associated with the \eqn{A} matrix of, \code{expmtA} is \eqn{exp(-A*(tree height))}, 
optimal regression is the \eqn{A^{-1}B}{\code{-solve(A)B}} matrix (if \eqn{A} is invertible, otherwise this will not exist),
\code{mPsi.rotated} is each of the regime effects multiplied by \eqn{1-\exp(-A*(tree height))}{1-exp(-A*(tree height))}, 
\code{cov.matrix} is the trait vector covariance matrix at the tree's height, 
\code{corr.matrix} is the trait vector correlation matrix at the tree's height, 
\code{conditional.cov.matrix} is the conditional covariance matrix of the OU type variables 
on the Brownian motion type at the tree's height, i.e. Cov[Y|X](tree height),
\code{conditional.corr.matrix} is the conditional correlation matrix of the OU type variables 
on the Brownian motion type at the tree's height, i.e. Corr[Y|X](tree height),
\code{stationary.cov.matrix} is the limit of the \code{conditional.cov.matrix}, \cr
\code{stationary.corr.matrix} is the limit of the \code{conditional.corr.matrix}, 
\code{optima.cov.matrix}
is the covariance matrix of the optimal process at the tree's height 
equalling \cr 
\eqn{(tree height)* A^{-1}B\Sigma_{xx}\Sigma_{xx}^{T}B^{T}A^{-T}}{\code{(tree height)* solve(A)\%*\%B\%*\%Sxx\%*\%t(Sxx)\%*\%t(B)\%*\%t(solve(A))}},
\code{optima.corr.matrix}
is the correlation matrix of the optimal process at time the tree's height,
\code{cov.with.optima} is the covariance matrix between the optimal process and the Y type variables process,
\code{corr.with.optima} is the correlation matrix between the optimal process and the Y type variables process, \cr
\code{evolutionary.regression} is the regression coefficient of E[Y|X](tree height).
Everything concerning the optimal process assumes A has positive real-part eigenvalues (in particular
it is invertible). Otherwise these will not exist.
\code{StS} is the infinitesimal covariance matrix,
\code{LogLik} the log--likelihood, dof the degrees of freedom, \code{m2loglik} is \eqn{-2}{-2}log--likelihood,
\code{aic} is the Akaike information criterion, \code{aic.c} is the Akaike information criterion corrected for small
sample size, sic is the Schwarz information criterion, \code{bic} is the Bayesian information criterion
(which is the same as the Schwarz information criterion) and \code{RSS} is the residual sum of squares.
The field \code{RSS_non_phylogenetic} is a residual sum of squares calculated without correcting
for the phylogeny--induced between species correlations, while the extension
conditional_on_predictors indicates that we consider the RSS for the variables
labelled as responses conditioned on the remaining variables. The \code{R2_phylaverage} field is
R2, where the alternative model is the phylogenetically weighted sample average 
(see \code{\link{OU_phylreg}}).
The last field \code{LogLik} is the log--likelihood at the point.

From version \code{2.0.0} of \pkg{mvSLOUCH} the data has to be passed as a matrix.
To underline this the data parameter's name has been changed to \code{mData}. 

From version \code{2.0.0} of \pkg{mvSLOUCH} the parameter \code{calcCI} has been removed. 
The package now offers the possibility of bootstrap confidence intervals, see
function \code{parametric.bootstrap}.
}
\value{
\item{FinalFound}{The point where the search procedure stopped. See Details for the description.}
\item{MaxLikFound}{The point with the highest likelihood found by the search procedure, 
if it is the same as the final point then this field equals "Same as final found". }
}
\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.
  
Butler, M.A. and A.A. King (2004) 
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.  
  American Naturalist 164:683-695.

Hansen, T.F. (1997) 
  Stabilizing selection and the comparative analysis of adaptation. 
  Evolution 51:1341-1351.

Hansen, T.F. and Bartoszek, K. (2012)
  Interpreting the evolutionary regression: the interplay between observational and 
  biological errors in phylogenetic comparative studies.
  Systematic Biology 61(3):413-425.

Hansen, T.F. and Pienaar, J. and Orzack, S.H. (2008) 
  A comparative method for studying adaptation to randomly evolving environment. 
  Evolution 62:1965-1977.

Labra, A., Pienaar, J. & Hansen, T.F. (2009) 
  Evolution of thermophysiology in Liolaemus lizards: adaptation, phylogenetic inertia and niche tracking. 
  The American Naturalist 174:204-220.

Mitov, V. and Bartoszek, K. and Asimomitis, G. and Stadler, T. (2020) 
  Fast likelihood calculation for multivariate Gaussian phylogenetic models with shifts
  Theoretical Population Biology 131:66-78.

Pienaar et al (in prep) An overview of comparative methods for
     testing adaptation to external environments.
}
\author{Krzysztof Bartoszek}
\note{
The slouch package is a recommended alternative if one has only a single response (Y) trait.
}

\section{Warning}{The estimation can take a long time and should be repeated a couple of times 
so that it is run from different starting positions.
The function can produce (a lot of) warnings and errors during the search procedure, 
this is nothing to worry about.
}

\seealso{
\code{\link[PCMBase]{PCMLik}}, \code{slouch::model.fit}, \code{\link{SummarizeMVSLOUCH}}, 
\code{\link{simulMVSLOUCHProcPhylTree}}, \cr \code{\link{parametric.bootstrap}}, 
\code{\link[stats]{optim}}
}

\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape. 
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(3)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

## 2 regimes
### Define a vector of regimes.
## regimes<-c("small","small","large","small")
## OUBMparameters<-list(vY0=matrix(1,ncol=1,nrow=1),A=matrix(0.5,ncol=1,nrow=1),
## B=matrix(2,ncol=1,nrow=1),mPsi=cbind("small"=1,"large"=-1),
## Syy=matrix(2,ncol=1,nrow=1),vX0=matrix(0,ncol=1,nrow=1),Sxx=diag(2,1,1),
## Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1))
## single regime for speed on CRAN
regimes<-c("small","small","small","small")
OUBMparameters<-list(vY0=matrix(1,ncol=1,nrow=1),A=matrix(0.5,ncol=1,nrow=1),
B=matrix(2,ncol=1,nrow=1),mPsi=cbind("small"=1),
Syy=matrix(2,ncol=1,nrow=1),vX0=matrix(0,ncol=1,nrow=1),Sxx=diag(2,1,1),
Syx=matrix(0,ncol=1,nrow=1),Sxy=matrix(0,ncol=1,nrow=1))


### Now simulate the data.
OUBMdata<-simulMVSLOUCHProcPhylTree(phyltree,OUBMparameters,regimes,NULL)
OUBMdata<-OUBMdata[phyltree$tip.label,,drop=FALSE]

### Try to recover the parameters of the mvOUBM model.
### maxiter here set to minimal working possibility, in reality it should be larger
### e.g. default of c(10,50,100)
### Also the Atype and Syytype variables should be changed, here set as simplest
### for speed of evaluation, e.g. Atype="DecomposablePositive", Syytype="UpperTri"
OUBMestim<-mvslouchModel(phyltree,OUBMdata,1,regimes,Atype="SingleValueDiagonal",
Syytype="SingleValueDiagonal",diagA="Positive",maxiter=c(1,2,1))
RNGversion(as.character(getRversion()))

\dontrun{ ##It takes too long to run this
## take a less trivial setup
phyltree<-ape::rtree(5)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

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
numboot=5,Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive")
}
}
\keyword{models}
\keyword{multivariate}
\keyword{htest}

