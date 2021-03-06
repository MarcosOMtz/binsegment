% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segtree.R
\name{fork}
\alias{fork}
\alias{fork.segtree}
\title{Fork Terminal Nodes}
\usage{
fork(x, ...)

\method{fork}{segtree}(tree, leaf, segvar, names, fast = tree$fast)
}
\arguments{
\item{tree}{An object of class \code{segtree}}

\item{leaf}{Either the name of a node or an actual node of \code{tree} which
is to be splitted}

\item{segvar}{Name of the segmentation variable to used to split}

\item{names}{A length-2 character vector containing *unique* names for the
resulting nodes}

\item{fast}{Whether (\code{FALSE}) or not (\code{TRUE}) to try out the splits
of the resulting nodes}
}
\value{
A tree with two new nodes corresponding to the new split that was requested
}
\description{
Actually performs a split on a given (usually terminal) node of a
\code{segree}.
}

