\name{generate.model.setups}
\alias{generate.model.setups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a list of model setups for the function
\code{\link{estimate.evolutionary.model}}.
}
\description{
The function generates a list of models that will be used by the 
function \cr \code{\link{estimate.evolutionary.model}}. A minimum
example list will be \code{list(list(evolmodel="bm"))}.
}
\usage{
generate.model.setups()
}
\details{
The function should really be a hidden one but is left available for the user as an
example how such a list of models should be generated.

The setting \code{Atype="Any"} means that one assumes the matrix \code{A} is eigendecomposable.
If \code{A} is defective, then the output will be erroneous.

None of the "signs" options for the model is generated, see the description of 
\code{mvslouchModel} and \code{ouchModel}.
}
\value{
A list with different models is returned. Each element of the list is a list with the following 
fields.
\itemize{
\item{\code{evolmodel} }{The evolutionary model, it may take one of the three values \code{"BM"} 
(Brownian motion model), \code{"ouch"} (OUOU model), \code{"mvslouch"} (OUBM model).}
\item{\code{Atype} }{The class of the \code{A} matrix, ignored if \code{evolmodel} equals 
\code{"BM"}. Otherwise it can take one of the following values: 
\code{"SingleValueDiagonal"}, \code{"Diagonal"},
\code{"UpperTri"}, \code{"LowerTri"}, \code{"SymmetricPositiveDefinite"},
\code{"Symmetric"}, \code{"DecomposablePositive"}, \cr \code{"DecomposableNegative"}, 
\code{"DecomposableReal"}, \code{"Invertible"},  \code{"TwoByTwo"}, \code{"Any"}.
}
\item{\code{Syytype} }{The class of the \code{A} matrix, ignored if \code{evolmodel} equals 
\code{"BM"}. Otherwise it can take one of the following values: 
\code{"SingleValueDiagonal"}, \code{"Diagonal"},
\code{"UpperTri"}, \code{"LowerTri"},  \code{"Symmetric"}, \code{"Any"}. 
}
\item{\code{diagA} }{Should the diagonal of \code{A} be forced to be positive (\code{TRUE}),
negative (\code{FALSE}) or the sign free to vary (\code{NULL})}
}
}

\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.
}
\author{Krzysztof Bartoszek}

\seealso{
\code{\link{estimate.evolutionary.model}}, \code{\link{mvslouchModel}}, \code{\link{ouchModel}} 
}
\examples{
model_setups<-generate.model.setups()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
