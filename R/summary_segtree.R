
#' Summarize Segmentation Trees
#'
#' Show the most important facts about a tree, including a text drawing of its
#' structure.
#'
#' @param tree tree An object of class \code{segtree}
#' @param s An object of class \code{summary.segtree}
#' @param details The level of details to show in the tree structure (see
#'   \code{\link{structure.segtree}})
#' @return An object of class \code{summary.segtree} which prints the most
#'   important features and structure of the tree in a nice and easy-to-see way
#' @export
summary.segtree <- function(tree){
  out <- list(
    population = nrow(tree$data),
    leaves = sum(sapply(tree$leaves, function(l) l$terminal)),
    intermediate_nodes = length(tree$leaves) - sum(sapply(tree$leaves, function(l) l$terminal)),
    formula = tree$formula,
    segvars = tree$segvars,
    gini = tree$gini,
    structure = structure.segtree(tree)
  )
  out$info <- data.frame(
    leaf = sapply(tree$leaves, function(l) l$name)
  )
  if(!tree$fast){
    out$info <- data.frame(
      leaf = out$info$leaf,
      depth = sapply(tree$leaves, function(l) length(l$segvars)),
      terminal = sapply(tree$leaves, function(l) ifelse(l$terminal, 'Yes', 'No')),
      population = sapply(tree$leaves, function(l) l$splits$population),
      p_population = sapply(tree$leaves, function(l) l$splits$population)/out$population,
      best_split_distr_A = sapply(tree$leaves, function(l) l$splits$table$p_pob_A[1]),
      best_split_distr_B = sapply(tree$leaves, function(l) l$splits$table$p_pob_B[1]),
      p_pos = sapply(tree$leaves, function(l) l$splits$p_pos_TOT),
      gini = sapply(tree$leaves, function(l) l$splits$gini_TOT),
      best_split_gini = sapply(tree$leaves, function(l) l$splits$table$gini_A_B[1])
    )
  }
  rownames(out$info) <- NULL
  class(out) <- 'summary.segtree'
  out
}

#' @rdname summary.segtree
#' @export
print.summary.segtree <- function(s, details = c(1, 2, 3, 0)){
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ >>\nSegmentation Tree Summary\n\nNumber of Leaves: %d\nNumber of Intermediate Nodes: %d\nTarget: %s\nRegression Variables:\n\t%s\nAvailable Segmentation Variables:\n\t%s\nPopulation: %s\nGlobal Gini Index: %.1f\n\nStructure:\n',
              s$leaves,
              s$intermediate_nodes,
              as.character(s$formula[2]),
              as.character(s$formula[3]),
              paste(s$segvars, collapse = ', '),
              format(s$population, scientific = F, big.mark = ','),
              100*s$gini))
  print.structure.segtree(s$structure, prefix='  ', details=details[1])
  cat('\n\nLeaf summary:\n\n')

  s$info$best_split_distr <- sprintf('(%.0f%% / %.0f%%)',
                                     100*s$info$best_split_distr_A,
                                     100*s$info$best_split_distr_B)
  s$info$population <- format(s$info$population, scientific = F, big.mark = ',')
  p_vars <- c('p_population','p_pos')
  gini_vars <- c('gini','best_split_gini')
  s$info[p_vars] <- apply(s$info[p_vars], 2, function(x) percent(x,1))
  s$info[gini_vars] <- 100*s$info[gini_vars]
  print(s$info[c('leaf','depth','terminal','population','p_population','best_split_distr',
                 'p_pos','gini','best_split_gini')], digits=3)
  cat('<< ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
}
