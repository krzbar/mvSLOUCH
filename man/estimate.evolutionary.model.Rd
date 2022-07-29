\name{estimate.evolutionary.model}
\alias{estimate.evolutionary.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper function to find best (out of BM, OU, OUOU, OUBM) fitting evolutionary model and 
estimate its parameters.
}
\description{
The \code{estimate.evolutionary.model} function calls the
\code{BrownianMotionModel},  \code{ouchModel} and \code{mvslouchModel}
functions with different classes of evolutionary model parameters.
It then compares the resulting estimates by the AICc (or BIC if AICc fails) and 
returns the best overall model. The user is recommended to install the suggested package 
\pkg{PCMBaseCpp} which significantly speeds up the calculations (see Details).
}
\usage{
estimate.evolutionary.model(phyltree, mData, regimes = NULL, 
root.regime = NULL, M.error = NULL, repeats = 5, model.setups = NULL, 
predictors = NULL, kY = NULL, doPrint = FALSE, pESS=NULL, 
estimate.root.state=FALSE, min_bl = 0.0003, maxiter=c(10,50,100))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phyltree}{
The phylogeny in \code{phylo} format. The tree can be obtained from e.g. a \code{nexus} file 
by the \code{read.nexus()} function from the \pkg{ape} package. The "standard" \pkg{ape} node
indexing is assumed: for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} 
and the root index \code{n+1}.
}
  \item{mData}{
A matrix with the rows corresponding to the tip species while the columns correspond to the traits.
The rows should be named by species \cr (field \code{phyltree$tip.label}), if not, then a warning
is thrown and the order of the species is assumed to be the same as the order in which the 
species are in the phylogeny (i.e. correspond to the node indices \code{1:n}, 
where \eqn{n}{n} is the number of tips). The columns should be named by traits, otherwise
a warning is thrown and generic names are generated.
}
  \item{regimes}{
A vector or list of regimes. If vector then each entry corresponds to each of \code{phyltree}'s
branches, i.e. to each row of \code{phyltree$edge}.If list then each list entry corresponds to a tip node
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
From version \code{2.0.0} of \pkg{mvSLOUCH} it is impossible to pass a single joint measurement error matrix 
for all the species and traits. 
}
  \item{repeats}{
How many starting points for the numerical maximum likelihood procedure should be tried for 
each model setup. On the first repeat for OUOU and OUBM modes the functions takes as the starting
point (for \code{A} and \code{Syy}) values based on the sample covariance matrix estimate, 
motivated by Bartoszek \eqn{\&}{&} Sagitov (2015)'s univariate results.
}

\item{model.setups}{
What models to try when searching for the best evolutionary model. 
This field may remain \code{NULL}, in this situation the function
generates using 

\code{.generate.basic.model.setups()} a basic list of models.
Allowed values are 
\itemize{
\item{\code{"basic"} }{A basic list of models to try out is generated, defined using 

\code{.generate.basic.model.setups()}. This list should be usually enough.}
\item{\code{"fundamental"} }{A slightly extended list of models to try out is generated, 
defined using \code{.generate.fund.model.setups()}. Compared to \code{"basic"} a few more 
models are added.}
\item{\code{"extended"} }{An extension of the \code{"fundamental"} list of models to try out. 
Defined using \code{.generate.ext.model.setups()} which at the moment calls \cr
\code{generate.model.setups()}.}
\item{\code{"all"} }{All possible models are generated, using 
\cr \code{.generate.all.model.setups()}. Running it will take an intolerable \cr amount of time.}
}
Alternatively the user is also free to provide their own list of models in this variable. 
Each element of the list is a list with fields.
\itemize{
\item{\code{evolmodel} }{The evolutionary model, it may take one of the three values \code{"bm"} 
(Brownian motion model), \code{"ouch"} (OUOU model), \code{"mvslouch"} \cr (OUBM model).}
\item{\code{Atype} }{The class of the \code{A} matrix, ignored if \code{evolmodel} equals 
\code{"bm"}. Otherwise it can take one of the following values: 
\code{"SingleValueDiagonal"}, \code{"Diagonal"},
\code{"UpperTri"}, \code{"LowerTri"}, \code{"SymmetricPositiveDefinite"},
\code{"Symmetric"}, \code{"DecomposablePositive"}, \code{"DecomposableNegative"}, \cr
\code{"DecomposableReal"}, \code{"Invertible"}, \code{"Any"}.
}
\item{\code{Syytype} }{The class of the \code{Syy} matrix, ignored if \code{evolmodel} equals 
\code{"bm"}. Otherwise it can take one of the following values: 
\code{"SingleValueDiagonal"}, \code{"Diagonal"},
\code{"UpperTri"}, \code{"LowerTri"},  \code{"Symmetric"}, \code{"Any"}. 
}
\item{\code{diagA} }{Should the diagonal of \code{A} be forced to be positive (\code{"Positive"}),
negative (\code{"Negative"}) or the sign free to vary (\code{NULL}). However, setting this to
a non-\code{NULL} value when \code{evolmodel} is \code{"mvslouch"} might be 
(but simulations concerning this are not conclusive) slightly detrimental to the optimization process 
if \code{Atype} is \code{"DecomposablePositive"}, \code{"DecomposableNegative"}, or
\code{"DecomposableReal"}. In these cases \code{A} is parametrized by its eigendecomposition.
Additional exponentiation of the diagonal, to ensure positivity, could (but this is uncertain)
make the exploration of the likelihood surface more difficult. The user is advised to
also try \code{diag=NULL}. In the case of \code{Atype} being \code{"SymmetricPositiveDefinite"},
the diagonal is always guaranteed to be positive.
}
\item{\code{signsA} }{WARNING: ONLY use this if you know what you are doing. Ignored if \code{evolmodel} 
equals \code{"bm"}. This allows the user to specify which elements of \code{A} are to be positive, negative 
or equal to specific values.See \code{\link{ouchModel}}  and \code{\link{mvslouchModel}} for a more specific 
description and important warnings.}
\item{\code{signsSyy} }{WARNING: ONLY use this if you know what you are doing. Ignored if \code{evolmodel} 
equals \code{"bm"}. This allows the user to specify which elements of \code{Syy} are to be positive, negative 
or equal to specific values. See \code{\link{ouchModel}}  and \code{\link{mvslouchModel}} for a more specific 
description and important warnings.}
\item{\code{signsB} }{WARNING: ONLY use this if you know what you are doing. Ignored if \code{evolmodel} 
does not equals \code{"mvslouch"}. This allows the user to specify which elements of \code{B} are to be 
positive, negative or equal to specific values. See \code{\link{mvslouchModel}} 
for a more specific description and important warnings.}
\item{\code{signsmPsi} }{WARNING: ONLY use this if you know what you are doing. Ignored if \code{evolmodel} 
equals \code{"bm"}. This allows the user to specify which elements of \code{mPsi} are to be positive, negative 
or equal to specific values.See \code{\link{ouchModel}} and \code{\link{mvslouchModel}} for a more specific 
description and important warnings.}
\item{\code{signsvY0} }{WARNING: ONLY use this if you know what you are doing. Ignored if \code{evolmodel} 
equals \code{"bm"}. This allows the user to specify which elements of \code{vY0} are to be positive, negative 
or equal to specific values.See \code{\link{ouchModel}}  and \code{\link{mvslouchModel}} for a more specific 
description and important warnings.}
\item{start_point_for_optim}{
A named list with starting parameters for of the parameters for be optimized by \code{optim()}, currently 
only \code{A} and \code{Syy} for \cr \code{evolmodel} equalling \code{"ouch"} or \code{"mvslouch"}. 
One may provide both or only one of them.
Make sure that the parameter is consistent with the other parameter restrictions as no check is done and this 
can result in undefined behaviour. For example one may provide this as (provided dimensions and other
parameter restrictions agree) \preformatted{start_point_for_optim=list(A=rbind(c(2,0),(0,4)), 
Syy=rbind(c(1,0.5),c(0,2))).}
 \item{parscale}{
A vector to calculate the \code{parscale} argument for \code{optim}. It is a named
vector with 3 entries, e.g. \cr \code{c("parscale_A"=3,"logparscale_A"=5,"logparscale_other"=1)}.
\cr The entry \code{parscale_A} is the scale for entries of the \code{A} matrix, \cr
\code{logparscale_A} is the scale for entries of the \code{A} matrix that are 
optimized over on the logarithmic scale, e.g. if eigenvalues are assumed to be positive,
then optimization is done over \code{log(eigenvalue)} for \code{A}'s eigendecomposition
and \code{logparscale_other} is the scale for entries other then of \code{A}
that are done on the logarithmic scale (e.g. \code{Syy}'s diagonal, or other entries
indicated as positive via \code{parameter_signs}). If not provided (or if a name of the vector
is misspelled), then made equal  to the example value provided above. For other elements, then 
mentioned above, that are optimized over by \code{optim()}, \code{1} is used for \code{optim()}'s 
\code{parscale}. It is advised that the user experiments with a couple of different values and 
reads \code{\link[stats]{optim}}'s man page. 
}
}
\item{estimateBmethod}{Only relevant for the OUBM models (optional), should \code{B} be estimated by maximum likelihood 
(default if not provided) value \code{"ML"} or generalized least squares (value \code{"GLS"}).
}

}
A minimum example list is \code{list(list(evolmodel="bm"))}. The functions that automatically generate
different types of models do NOT use any of the "signs" parameters. Hence, in these models 
all parameters (under the appropriate parametrization) will be free to vary.
}
  \item{predictors}{
A vector giving the numbers of the columns from
\code{dfdata} which are to be considered predictor ones, \emph{i.e.} conditioned
on in the program output. 
A vector giving the numbers of the columns from \code{mData} matrix which are to be considered predictor ones, 
\emph{i.e.} conditioned on in the program output. If not provided then in for the OUBM model the columns \code{(kY+1):ncol(mData)}, i.e. the 
"BM" ones, are treated as predictors. Otherwise, none will be considered to be predictors.
}
  \item{kY}{
Number of "Y" (response) variables, for the OUBM models.
The first \code{kY}
columns of \code{mY} are the "OU" ones, while the rest the "BM" ones.
In more detail this value determines the number of columns of the \code{mData} matrix to treat as response variables ("OU" ones). For example, a 
value of 1 means that only the first column is treated as a response variable, while a value of 3 means the first three 
columns are treated as response variables. Any predictor variables ("BM" ones) the user is interested in setting for a particular model 
should therefore be placed in the final columns of the \code{mData} matrix, allowing for selecting select \code{kY} columns before 
this as response variables ("OU" ones).
}
  \item{doPrint}{
Should the function print out information on what it is doing (TRUE) or keep silent 
(default FALSE).
}
  \item{pESS}{
Should the function also find the best model taking into account the phylogenetic effective sample size
and it so what method. If \code{NULL}, then do not take this into account. Otherwise one
of \code{"reg"} ("regression" effective sample size that takes into account all of the 
correlations between species explicitly),
\code{"mean"} (mean value effective sample size \eqn{1^{T}R^{-1}1}{\code{t(1)\%*\%solve(R)\%*\%1}}, where \eqn{R} 
is the interspecies correlation matrix),
\code{"MI"} (mutual information effective sample size), 
\code{"mvreg"}(multivariate version of "regression" effective sample size when each species is 
described by a suite of traits) , 
\code{"mvMI"} (multivariate mutual information effective sample size when each species is described by 
a suite of traits) indicating
the way to calculate the pESS. The default (\code{NULL}) is not to do any pESS calculations as these will be slow. 
They require the construction of the between-species-between-traits variance covariance matrix and
hence do not fully take advantage of the speed-up offered by \pkg{PCMBase}. If \code{pESS="only_calculate"},
then all possible pESS values are calculated but no model selection is done based on them.
}
\item{estimate.root.state}{
    Should the root state be estimate \code{TRUE} (not recommended) or set at the optimum 
    \code{FALSE} (recommended). Root state estimation is usually unreliable hence
    if fossil measurements are available prediction based on them and the estimated
    model will probably be more accurate. If there is only one regime, then 
    estimation of the root state separately is impossible and will not be allowed.
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
  in the Brownian motion and measurement error case.
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

If \code{model.setups} is left at the default value the 
function will take a long time to run, as it performs estimation for each model 
(\code{generate.model.setups} generates 90 setups) times the value in repeats. 
Therefore if the user has particular hypotheses in mind then it is advisable to
prepare their own list.
If the \code{Syy} matrix is assumed to be upper-triangular and the starting conditions
based on Bartoszek \eqn{\&}{&} Sagitov (2015)'s results are used then the factorization
of \eqn{\Sigma=\Sigma_{yy} \Sigma_{yy}^{T}}{\code{S=Syy\%*\%t(Syy)}} into \eqn{\Sigma_{yy}}{Syy} is 
done using the procedure described in 
\url{https://math.stackexchange.com/questions/2039477/cholesky-decompostion-upper-triangular-or-lower-triangular}.

From version \code{2.0.0} of \pkg{mvSLOUCH} the data has to be passed as a matrix.
To underline this the data parameter's name has been changed to \code{mData}. 

If AICc fails, then the function will use BIC to select between models. This is extremely 
unlikely essentially only when AICc is infinite, i.e. the model is saturated 
(number of observations equals number of data points).
}
\value{
A list is returned that describes the results of the search. See the help for
\code{\link{BrownianMotionModel}}, \code{\link{ouchModel}} and
\code{\link{mvslouchModel}} for the description of the lower level entries.
The elements of this list are the following
\item{\code{BestModel} }{The resulting best model found. Included are the model parameters, 
a "first-glance" qualitative description of the model, the most important parameters of the 
process (half-lives and regressions in the case of OU models) and what to call to obtain 
standard errors. It takes a long time to obtain them so calculating them is not part of the 
standard procedure.}
\item{\code{BestModelESS} }{Only if \code{pESS} was \code{TRUE}. The resulting best model found
taking into account the phylogenetic essential sample size. Included are the model parameters, 
a "first-glance" qualitative description of the model, the most important parameters of the 
process (half-lives and regressions in the case of OU models) and what to call to obtain 
bootstrap confidence intervals. It takes a long time to obtain them so calculating them is not part of the 
standard procedure.}
\item{\code{testedModels} }{A list of results for each tried model.}
\item{\code{model.setups} }{A list of models tried.}
\item{\code{repeats} }{How many starting points were tried per model.} 
}

\references{
Ane, C. (2008) 
  Analysis of comparative data with hierarchical autocorrelation. 
  Annals of Applied Statistics 2:1078-1102.

Bartoszek, K. (2016)
  Phylogenetic effective sample size.
  Journal of Theoretical Biology 407:371-386.

Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.

Bartoszek, K. and Sagitov, S. (2015) 
  Phylogenetic confidence intervals for the optimal trait value. 
  Journal of Applied Probability 52(4):1115-1132.
  
Butler, M.A. and A.A. King (2004) 
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution.  
  American Naturalist 164:683-695.

Hansen, T.F. and Pienaar, J. and Orzack, S.H. (2008) 
  A comparative method for studying adaptation to randomly evolving environment. 
  Evolution 62:1965-1977.

Mitov, V. and Bartoszek, K. and Asimomitis, G. and Stadler, T. (2020) 
  Fast likelihood calculation for multivariate Gaussian phylogenetic models with shifts
  Theoretical Population Biology 131:66-78.

Xiao, H and Bartoszek, K. and Lio P. (2018)
  Multi--omic analysis of signalling factors in inflammatory comorbidities. 
  BMC Bioinformatics, Proceedings from the 12th International BBCC conference 
  19:439.  
}
\author{Krzysztof Bartoszek}

\note{
The engine behind the likelihood calculations is called from \pkg{PCMBase}.
The \code{slouch} package is a recommended alternative if one has a OUBM models and 
only a single response (Y) trait. 
The \code{mvMORPH}, \code{ouch} and \code{Rphylpars} packages consider multivariate OU models 
and looking at them could be helpful.
}

\seealso{
\code{\link[ouch]{brown}}, \code{\link[mvMORPH]{mvBM}}\code{\link{BrownianMotionModel}}, \code{\link{SummarizeBM}}, 
\code{\link{simulBMProcPhylTree}}, \code{\link[ouch]{hansen}}, \code{\link[mvMORPH]{mvOU}}, \cr \code{\link{ouchModel}}, 
\code{\link{SummarizeOUCH}}, \code{\link{simulOUCHProcPhylTree}},
\code{slouch::model.fit}, \code{\link[PCMBase]{PCMLik}}, \cr \code{\link{mvslouchModel}},
\code{\link{SummarizeMVSLOUCH}}, \code{\link{simulMVSLOUCHProcPhylTree}},\cr  \code{\link{parametric.bootstrap}},
\code{\link[stats]{optim}}
}
\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
### We will first simulate a small phylogenetic tree using functions from ape.
### For simulating the tree one could also use alternative functions, e.g. sim.bd.taxa 
### from the TreeSim package
phyltree<-ape::rtree(4)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

### Define a vector of regimes.
regimes<-c("small","small","large","small","large","small")

### Define SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),"large"=c(-1,1,0.5)),
Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

### Now simulate the data.
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

## set up for a trivial, single model setup case (for running time)
## in a real analysis you should carefully choose between what models
## you want to do model selection
model_setups<-list(list(evolmodel="bm"))

### Try to recover the parameters of the OUOU model.
### maxiter here set to minimal working possibility, in reality it should be larger
### e.g. default of c(10,50,100)
estimResults<-estimate.evolutionary.model(phyltree,OUOUdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=1,model.setups=model_setups,predictors=c(3),
kY=2,doPrint=TRUE,pESS=NULL,maxiter=c(1,1,1))

### After this step you can look at the best estimated model and use the
### parametric.bootstrap() function to obtain bootstrap confidence intervals
RNGversion(as.character(getRversion()))

\dontrun{ ##It takes too long to run this
## take a less trivial setup
phyltree<-ape::rtree(5)

## The line below is not necessary but advisable for speed
phyltree<-phyltree_paths(phyltree)

### Define a vector of regimes.
regimes<-c("small","small","large","small","small","large","large","large")

### Define SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(c(1,-1,0.5),nrow=3,ncol=1),
A=rbind(c(9,0,0),c(0,5,0),c(0,0,1)),mPsi=cbind("small"=c(1,-1,0.5),"large"=c(-1,1,0.5)),
Syy=rbind(c(1,0.25,0.3),c(0,1,0.2),c(0,0,1)))

### Now simulate the data.
OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[phyltree$tip.label,,drop=FALSE]

## set up for two very simple (for example usage) models to compare between
## in a real analysis you should carefully choose between what models
## you want to do model selection, the default 
## model_setups<-NULL provides a wide selection of models

model_setups<-list(list(evolmodel="bm"),list(evolmodel="ouch",
"Atype"="SingleValueDiagonal","Syytype"="SingleValueDiagonal","diagA"="Positive"))

### Try to recover the parameters of the OUOU model.
estimResults<-estimate.evolutionary.model(phyltree,OUOUdata,regimes=regimes,
root.regime="small",M.error=NULL,repeats=3,model.setups=model_setups,predictors=c(3),
kY=2,doPrint=TRUE,pESS=NULL,maxiter=c(10,50,100))
}
}
\keyword{models}
\keyword{multivariate}
\keyword{htest}
