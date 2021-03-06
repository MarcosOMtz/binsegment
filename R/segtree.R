
require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)

#' Create Segmentation Trees
#'
#' \code{segtree} builds a basic segmentation tree whose nodes can then be splitted with \code{fork} and evaluated with a variety of functions.
#'
#' @param formula A formula specifying the response and features for the logistic regressions
#' @param data A \code{data.frame} containing the variables specified in formula and in \code{segvars}
#' @param segvars A character vector with the names of the segmentation variables to be tried
#' @param fast If \code{FALSE} (the default), evaluate the splits for each node. Otherwise construct a simpler tree which takes less time to build
#' @return An object of class \code{segtree}
#' @seealso \code{\link{fork}}, \code{\link{summary.segtree}}, \code{\link{performance.segtree}}, \code{\link{split_data}}, \code{\link{split.segtree}}
#' @export
segtree <- function(formula, data, segvars, fast=FALSE, ...){

  resp <- as.character(formula)[2]
  regvars <- strsplit(as.character(formula)[3], split = ' \\+ ')[[1]]
  if(!(resp %in% names(data))){
    stop(sprintf('The response variable %s is not in the dataset.', resp))
  }
  if(!all(regvars %in% names(data))){
    stop(sprintf('The following regression variables are not in the dataset: %s',
                 paste(regvars[!(regvars %in% names(data))], collapse=', ')))
  }
  if(!all(segvars %in% names(data))){
    stop(sprintf('The following segmentation variables are not in the dataset: %s',
                 paste(segvars[!(segvars %in% names(data))], collapse=', ')))
  }
  if(fast){
    splits <- list()
    g0 <- NULL
  } else{
    splits <- split_.formula(formula, data, segvars, ...)
    m0 <- glm(formula, data, family=binomial(link='logit'))
    yhat0 <- predict(m0, data, type='link')
    g0 <- optimbucket::performance(yhat0, m0$y, ...)$gini
  }
  out <- list(
    formula = formula,
    leaves = list('root'=leaf(segvars=NULL, levels=NULL, name='root',
                              splits=splits, terminal=TRUE)),
    data = data,
    segvars = segvars,
    gini = g0,
    fast = fast
  )
  class(out) <- c('segtree')
  out
}

#' @rdname segtree
#' @export
print.segtree <- function(tree, ...){
  nnodes <- length(tree$leaves)
  nleaves <- sum(sapply(tree$leaves, function(l) l$terminal))
  target <- as.character(tree$formula[2])
  feats <- as.character(tree$formula[3])
  segvars <- paste(tree$segvars, collapse = ', ')
  population <- format(nrow(tree$d), scientific = F, big.mark = ',')
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ >>\nSegmentation Tree\n\nNumber of Leaves: %d\nNumber of Intermediate Nodes: %d\nTarget: %s\nRegression Variables:\n\t%s\nAvailable Segmentation Variables:\n\t%s\nPopulation: %s\nGlobal Gini Index: %.1f\nLeaves:\n\n',
              nleaves,
              nnodes - nleaves,
              target,
              feats,
              segvars,
              population,
              100*tree$gini))
  for(l in tree$leaves){
    print(l, global_pop = nrow(tree$data))
  }
  cat('<< ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
}

#' Evaluate Splits on a Tree
#'
#' This function is used internally to try out all the remaining segmentation
#' variables on a given node.
#'
#' @param tree An object of class \code{segtree}
#' @param leaf Either the name of a node or an actual node of \code{tree}
#' @return A simple object of class \code{segtree.split} containing the results
#'   of splitting with each remaining segmentation variable
#' @export
split.segtree <- function(tree, leaf){
  if(class(tree) != 'segtree' ||
     (class(leaf) != 'leaf' &&
      class(leaf) != 'character')){
    stop('Please provide a segtree and a leaf.')
  } else if(is.character(leaf)){
    leaf_ix <- which(names(tree$leaves) == leaf)
    if(length(leaf_ix) == 0){
      stop(sprintf('Leaf %s does not exist.', leaf))
    } else{
      leaf <- tree$leaves[[leaf_ix]]
    }
  }

  if(leaf$name == 'root'){
    dat <- tree$data
    segvars <- tree$segvars
  } else{
    levs <- matrix(rep(leaf$levels, nrow(tree$data)),
                   nrow = nrow(tree$data),
                   ncol=length(leaf$levels),
                   byrow=T)
    ix <- (tree$data[,leaf$segvars] == levs) %>%
      apply(1, all)
    rm(levs)
    dat <- tree$data[ix,]
    segvars <- tree$segvars[!(tree$segvars %in% leaf$segvars)]
  }

  split_.formula(tree$formula, dat, segvars)
}

#' @rdname split.segtree
#' @export
print.segtree.split <- function(s, global_pop = NULL){
  if(is.null(global_pop)){
    cat(sprintf('Population: %s\nTotal P(Y = 1): %.2f%%\nTotal Gini Index: %.1f\nDetails:\n\n',
                format(s$population, big.mark = ',', scientific = F),
                100*s$p_pos_TOT, 100*s$gini_TOT))
    p_vars <- c('p_pob_A','p_pob_B','p_pos_A','p_pos_B')
    gini_vars <- c('gini_A_B', 'gini_A', 'gini_B')
    s$table[p_vars] <- apply(s$table[p_vars], 2, scales::percent)
    s$table[gini_vars] <- 100*s$table[gini_vars]
    print(s$table, digits = 3)
  } else{
    p_glob <- s$population/global_pop
    cat(sprintf('\nGlobal Population: %s\nPopulation: %s (%.1f%%)\nTotal P(Y = 1): %.2f%%\nTotal Gini Index: %.1f\nDetails:\n',
                format(global_pop, big.mark = ',', scientific = F),
                format(s$population, big.mark = ',', scientific = F),
                100*p_glob,
                100*s$p_pos_TOT,
                100*s$gini_TOT))
    s$table$p_glob_A <- s$table$p_pob_A*p_glob
    s$table$p_glob_B <- s$table$p_pob_B*p_glob
    p_vars <- c('p_pob_A','p_pob_B','p_pos_A','p_pos_B','p_glob_A','p_glob_B')
    gini_vars <- c('gini_A_B', 'gini_A', 'gini_B')
    s$table[p_vars] <- apply(s$table[p_vars], 2, function(x) percent(x,1))
    s$table[gini_vars] <- 100*s$table[gini_vars]
    print(s$table, digits = 3)
  }
}

#' Fork Terminal Nodes
#'
#' Actually performs a split on a given (usually terminal) node of a
#' \code{segree}.
#'
#' @param tree An object of class \code{segtree}
#' @param leaf Either the name of a node or an actual node of \code{tree} which
#'   is to be splitted
#' @param segvar Name of the segmentation variable to used to split
#' @param names A length-2 character vector containing *unique* names for the
#'   resulting nodes
#' @param fast Whether (\code{FALSE}) or not (\code{TRUE}) to try out the splits
#'   of the resulting nodes
#' @return A tree with two new nodes corresponding to the new split that was requested
#' @export
fork <- function(x, ...) UseMethod('fork')

#' @rdname fork
#' @export
fork.segtree <- function(tree, leaf, segvar, names=NULL, fast=tree$fast){
  if(is.character(leaf)){
    leaf_ix <- which(names(tree$leaves) == leaf)
    if(length(leaf_ix) == 0){
      stop(sprintf('Leaf %s does not exist.', leaf))
    } else{
      leaf <- tree$leaves[[leaf_ix]]
    }
  } else if(class(leaf) == 'leaf'){
    leaf_ix <- which(names(tree$leaves) == leaf$name)
    if(length(leaf_ix) == 0){
      stop(sprintf('Leaf %s does not exist.', leaf$name))
    }
  } else{
    stop('Please supply a leaf name or a leaf object.')
  }

  if(!is.character(segvar) |
     length(segvar) != 1 |
     !(segvar %in% tree$segvars)){
    stop(sprintf('Invalid value for segvar (%s). Please choose a variable in tree$segvars.', segvar))
  }

  if(!(segvar %in% tree$segvars)){
    warning(sprintf('Segmentation variable %s is not present in tree$segvars and so will not appear in any splits. Run a code similar to "<your_tree>$segvars <- c(<your_tree>$segvars, %s)" to add it and avoid this warning.', segvar, segvar))
  }
  classes <- levels(as.factor(tree$data[[segvar]]))[1:2]

  if(is.null(names) | !is.character(names) | length(names) < 2){
    if(is.null(names)){
      message('No names for new nodes provided. Choosing names automatically...')
    } else{
      warning('Improper names chosen. Choosing automatically...')
    }
    names <- paste(sprintf('%s > %s =', leaf$name, segvar), classes)
  }

  leaf_A <- leaf(
    segvars = c(leaf$segvars, segvar),
    levels = c(leaf$levels, classes[1]),
    name = names[1],
    terminal = TRUE
  )
  leaf_B <- leaf(
    segvars = c(leaf$segvars, segvar),
    levels = c(leaf$levels, classes[2]),
    name = names[2],
    terminal = TRUE
  )
  tree$leaves[[leaf_ix]]$children <- c(tree$leaves$children, names)
  tree$leaves[[leaf_ix]]$terminal <- FALSE
  new_leaves <- list(leaf_A,leaf_B)
  names(new_leaves) <- names[1:2]
  tree$leaves <- c(tree$leaves, new_leaves)
  if(!fast){
    tree$leaves[[names[1]]]$splits <- split.segtree(tree, names[1])
    tree$leaves[[names[2]]]$splits <- split.segtree(tree, names[2])
  }
  tree
}

# plot.segtree <- function(tree){
#   defs <- sapply(tree$leaves, function(ll){
#     def <- paste(paste(ll$segvars, ll$levels, sep=' = '), collapse=' > ')
#     cat(sprintf(''))
#   })
# }
# plot.segtree(tt2)


#' Evaluate the Performance of a Segmentation
#'
#' Compares the Gini Index before and after applying the segmentation.
#'
#' @param tree An object of class \code{segtree}
#' @param details The level of details to show in the tree structure (see
#'   \code{\link{structure.segtree}})
#' @return An object of class \code{performance.segtree} with the most important
#'   stats of the tree's performance
#' @export
performance <- function(x, ...) UseMethod('performance')

#' @rdname performance
#' @export
performance.segtree <- function(tree, ...){
  leaves <- tree$leaves[sapply(tree$leaves, function(l) l$terminal)]
  yhats <- list()
  ys <- list()
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\nCalculating Performance of Tree with %d Terminal Nodes...\n\n',
              length(leaves)))
  for(i in 1:length(leaves)){
    cat(sprintf('(%d/%d) Running Regression for Terminal Node %s\t@ %s\n',
                i, length(leaves), leaves[[i]]$name, Sys.time()))
    if(is.null(leaves[[i]]$levels)){
      dat <- tree$data
    } else{
      levs <- matrix(rep(leaves[[i]]$levels, nrow(tree$data)),
                     nrow = nrow(tree$data),
                     ncol=length(leaves[[i]]$levels),
                     byrow=T)
      ix <- (tree$data[,leaves[[i]]$segvars] == levs) %>%
        apply(1, all)
      rm(levs)
      dat <- tree$data[ix,]
    }
    m <- glm(tree$formula, data=dat, family=binomial(link='logit'))
    ys[[i]] <- m$y
    yhats[[i]] <- predict(m, newdata=dat, type='link')
  }
  y <- unlist(ys)
  yhat <- unlist(yhats)
  out <- list(
    tree = tree,
    performance = optimbucket::performance(yhat, y, ...)
  )
  class(out) <- 'performance.segtree'
  out
}

#' @rdname performance
#' @export
print.performance.segtree <- function(object, details = c(1, 2, 3, 0), ...){
  tree <- object$tree
  perf <- object$performance
  nnodes <- length(tree$leaves)
  nleaves <- sum(sapply(tree$leaves, function(l) l$terminal))
  target <- as.character(tree$formula[2])
  feats <- as.character(tree$formula[3])
  segvars <- paste(tree$segvars, collapse = ', ')
  population <- format(nrow(tree$d), scientific = F, big.mark = ',')
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ >>\nSegmentation Tree\n\nNumber of Leaves: %d\nNumber of Intermediate Nodes: %d\nTarget: %s\nRegression Variables:\n\t%s\nAvailable Segmentation Variables:\n\t%s\nPopulation: %s\nGlobal Gini Index: %.1f\n\n>>>>>>> Gini Index Using Segmentation: %.1f <<<<<<<\n\nStructure:\n',
              nleaves,
              nnodes - nleaves,
              target,
              feats,
              segvars,
              population,
              100*tree$gini,
              100*perf$gini))
  print(structure.segtree(tree, 'root'), details=details[1])
  cat('\n<< ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~\n')
}








