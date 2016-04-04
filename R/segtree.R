
require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)

segtree <- function(formula, data, segvars, fast=FALSE, ...){
  if(fast){
    splits <- list()
    g0 <- NULL
  } else{
    splits <- split_.formula(formula, data, segvars, ...)
    m0 <- glm(formula, data, family=binomial(link='logit'))
    yhat0 <- predict(m0, data, type='link')
    g0 <- performance(yhat0, m0$y, ...)$gini
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

print.segtree <- function(tree, ...){
  nnodes <- length(tree$leaves)
  nleaves <- sum(sapply(tree$leaves, function(l) l$terminal))
  target <- as.character(tree$formula[2])
  feats <- as.character(tree$formula[3])
  segvars <- paste(tree$segvars, collapse = ', ')
  population <- format(nrow(tree$d), scientific = F, big.mark = ',')
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ >>\nSegmentation Tree\n\nNumber of Leaves: %d\nNumber of Intermediate Nodes: %d\nTarget: %s\nRegression Variables:\n\t%s\nAvailable Segmentation Variables:\n\t%s\nPopulation: %s\nGlobal Gini Index: %.3f\nLeaves:\n\n',
              nleaves,
              nnodes - nleaves,
              target,
              feats,
              segvars,
              population,
              tree$gini))
  for(l in tree$leaves){
    print(l, global_pop = nrow(tree$data))
  }
  cat('<< ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
}

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

print.segtree.split <- function(s, global_pop = NULL){
  if(is.null(global_pop)){
    cat(sprintf('Population: %s\nTotal P(Y = 1): %.2f%%\nTotal Gini Index: %.3f\nDetails:\n\n',
                format(s$population, big.mark = ',', scientific = F),
                100*s$p_pos_TOT, s$gini_TOT))
    p_vars <- c('p_pob_A','p_pob_B','p_pos_A','p_pos_B')
    s$table[p_vars] <- apply(s$table[p_vars], 2, scales::percent)
    print(s$table, digits = 3)
  } else{
    p_glob <- s$population/global_pop
    cat(sprintf('\nGlobal Population: %s\nPopulation: %s (%.1f%%)\nTotal P(Y = 1): %.2f%%\nTotal Gini Index: %.3f\nDetails:\n',
                format(global_pop, big.mark = ',', scientific = F),
                format(s$population, big.mark = ',', scientific = F),
                100*p_glob,
                100*s$p_pos_TOT, s$gini_TOT))
    s$table$p_glob_A <- s$table$p_pob_A*p_glob
    s$table$p_glob_B <- s$table$p_pob_B*p_glob
    p_vars <- c('p_pob_A','p_pob_B','p_pos_A','p_pos_B','p_glob_A','p_glob_B')
    s$table[p_vars] <- apply(s$table[p_vars], 2, function(x) percent(x,1))
    print(s$table, digits = 3)
  }
}

fork.segtree <- function(tree, leaf, segvar, names, fast=tree$fast){
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

  classes <- levels(as.factor(tree$data[[segvar]]))[1:2]
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












