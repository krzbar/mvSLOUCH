\name{phyltree_paths}
\alias{phyltree_paths}
\title{Extract path information from a phylogenetic tree}
\description{
The function computes for each node its path to the root and
its distance to the root. It returns an ``enhanced'' 
\code{phylo} type tree.}

\usage{
phyltree_paths(phyltree)
}

\arguments{
  \item{phyltree}{The phylogeny - an object of class \code{phylo}, i.e. tree in \pkg{ape} format.
  The "standard" \pkg{ape} node indexing is assumed: for a tree with \eqn{n}{n} tips, 
  the tips should have indices \code{1:n} and the root index \code{n+1}.
  }
}

\details{The function removes a root edge, i.e. \code{$root.edge}
if one is present.
}

\value{The function returns a \code{phylo} type tree with the below additional
fields.
\item{Ntips}{Number of tips on the tree.}
\item{path.from.root}{A list of length equalling the number of nodes. Each 
entry is a list made up of two fields \code{nodes} and \code{edges}. 
\code{nodes} are the nodes on the path to the root and \code{edges} the edges.}
\item{time.of.nodes}{A vector of length equalling the number of nodes. Each 
entry is the node's distance from the root. This is only calculated if the input
tree has the \code{$edge.length} field.}
\item{tree_height}{The height of the tree if it is ultrametric, otherwise the 
length of the longest path from root to tip.}
\item{tip_species_index}{The node numbers corresponding to tip nodes, should equal \code{1:n}.}
\item{internal_nodes_index}{The node numbers corresponding to internal nodes, should equal
\code{(n+1):(2n-1)}.}
\item{root_index}{The node number corresponding to the root, should equal \code{n+1}.}
}

\author{Krzysztof Bartoszek}

\note{The \code{ape} and \code{phangorn} packages include related tree manipulation functions.}

\seealso{\code{\link{ape}}, \code{phangorn}}

\examples{
    RNGversion(min(as.character(getRversion()),"3.6.1"))
    set.seed(12345, kind = "Mersenne-Twister", normal.kind = "Inversion")
    phyltree<-ape::rtree(5)
    phyltree_augmented<-phyltree_paths(phyltree)
    RNGversion(as.character(getRversion()))
}

\keyword{manip}
