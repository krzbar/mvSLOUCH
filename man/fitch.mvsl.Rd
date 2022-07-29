% --- Source file: fitch.mvsl.R ---
\name{fitch.mvsl}
\alias{fitch.mvsl}
\title{Unordered Fitch parsimony reconstruction of discrete character states}
\description{
Implements an unordered Fitch parsimony reconstruction of discrete niche variables 
for use in the OU models where optima are modeled on discrete, categorical niche encodings. 
Allows for delayed and accelerated transformations to deal with ambiguities. 
Function was originally the \code{fitch()} function from the \pkg{slouch} package.
}
\usage{
fitch.mvsl(phyltree, niche, deltran = FALSE, acctran = FALSE, root = NULL)
}
\arguments{
  \item{phyltree}{The phylogenetic tree in \pkg{ape} (\code{phylo}).
  For a phylogeny in \code{phylo} format the "standard" \pkg{ape} node indexing is assumed: 
  for a tree with \eqn{n}{n} tips, the tips should have indices \code{1:n} and the 
  root index \code{n+1}.
  }
  \item{niche}{The specific niche variable in the \pkg{slouch} data.frame to be 
  reconstructed, entered as data.frame\code{$}niche . The order of the niche's regimes
  has to correspond to the order of the tip nodes in \code{phyltree}. }
  \item{deltran}{Implements a delayed transformation algorithm in order to 
  deal with ambiguous nodes}
  \item{acctran}{Implements an accelerated transformation algorithm to deal 
  with ambiguous nodes}
  \item{root}{An optional argument allowing the user to define a character 
  state for the root (useful if the root node is ambiguously reconstructed)}
}
\details{
The \code{fitch.mvsl} function is meant to be interactive, where the user acts on the advice 
given in the returned messages whilst attempting to reconstruct ancestral states. 
If the root node is ambiguous after an initial reconstruction (a message will be printed to 
the screen if this is the case), this needs to be set by the user using the \code{root = "state"} 
argument in the function call. Any remaining ambiguous nodes can then be dealt with by specifying 
\code{deltran} or \code{acctran ="TRUE"} in the function call 
}
\value{
The \code{fitch.mvsl} function returns a list with two or three elements. The first, 
\code{$branch_regimes} is a vector of reconstructed character states. Each entry of
the vector corresponds to the respective edge in the \code{$edge} field in the provided tree.
Notice that entries correspond to edges and not to nodes. If you require correspondence
with nodes, then you can treat the given edge entry as the value for the node ending the edge.
Actually, this is what the algorithm in the function estimates. 
The second field of the output object, \code{$root_regime}s is the regime at the root of the tree.
If the provided tree was a raw \code{phylo} object, then the function will also return 
an enhanced version of it (field \code{$phyltree}). This is the tree that results from calling
\pkg{mvSLOUCH}\code{::phyltree_paths(phyltree)} on the originally provided tree. This
enhanced version is returns as calculating it is costly and the user might want to re-use
it in some downstream analysis with \pkg{mvSLOUCH}. All \pkg{mvSLOUCH} user-level functions
first enhance the provided phylogeny by \pkg{mvSLOUCH}\code{::phyltree_paths()}, but they
first check if it is not already enhanced. 
}
\references{ 
Fitch, M.W. (1971) Defining the course of Evolution: Minimum change for a specific tree topology. 
\emph{Systematic Zoology} \bold{20}:406--416.

Swofford, D. L. and W.P. Maddison (1987) Reconstructing ancestral character states under 
Wagner parsimony. \emph{Mathematical Biosciences} \bold{87}: 199--229. 
}
\author{Jason Pienaar \email{jasonpienaar@gmail.com}}
\seealso{
\code{slouch::fitch}, \code{slouch::slouchtree.plot}, \code{slouch::model.fit}, 
\code{slouch::ouch2slouch}
}

\examples{
RNGversion(min(as.character(getRversion()),"3.6.1"))
set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
phyltree<-ape::rtree(5)

regimes<-c("A","B","B","C","C")
regimesFitch<-fitch.mvsl(phyltree,regimes,root=1,deltran=TRUE)
RNGversion(as.character(getRversion()))
}
\keyword{models}
\keyword{htest}
