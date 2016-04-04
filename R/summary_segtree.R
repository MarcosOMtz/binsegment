
summary.segtree <- function(tree){
  out <- list(
    population = nrow(tree$data),
    leaves = length(tree$leaves),
    formula = tree$formula,
    segvars = tree$segvars,
    gini = tree$gini
  )
  out$info <- data.frame(
    leaf = sapply(tree$leaves, function(l) l$name)
  )
  if(!tree$fast){
    out$info <- data.frame(
      leaf = out$info$leaf,
      depth = sapply(tree$leaves, function(l) length(l$segvars)),
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

print.summary.segtree <- function(s, ...){
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ >>\nSegmentation Tree Summary\n\nNumber of leaves: %d\nTarget: %s\nRegression Variables:\n\t%s\nAvailable segmentation variables:\n\t%s\nPopulation: %s\nGlobal Gini Index: %.3f\nLeaf summary:\n\n',
              s$leaves,
              as.character(s$formula[2]),
              as.character(s$formula[3]),
              paste(s$segvars, collapse = ', '),
              format(s$population, scientific = F, big.mark = ','),
              s$gini))
  s$info$best_split_distr <- sprintf('(%.0f%% / %.0f%%)',
                                     100*s$info$best_split_distr_A,
                                     100*s$info$best_split_distr_B)
  s$info$population <- format(s$info$population, scientific = F, big.mark = ',')
  p_vars <- c('p_population','p_pos')
  s$info[p_vars] <- apply(s$info[p_vars], 2, scales::percent)
  print(s$info[c('leaf','depth','population','p_population','best_split_distr',
                 'p_pos','gini','best_split_gini')], digits=3)
  cat('<< ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
}
