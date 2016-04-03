
require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)

split_one <- function(formula, data, segvar, ngroups = 100, ...){
  m0 <- glm(formula, data, family=binomial(link='logit'))
  yhat0 <- predict(m0, data, type='link')
  g0 <- performance(yhat0, m0$y, ngroups = ngroups, ...)$gini

  classes <- levels(as.factor(data[[segvar]]))[1:2]
  ix_A <- data[[segvar]] == classes[1]
  ix_B <- data[[segvar]] == classes[2]
  if(sum(ix_A) == 0 || sum(ix_B) == 0){
    warning(sprintf('Variable %s has no levels of some of its classes. Skipping...', segvar))
    return(NULL)
  }
  data_A <- data[ix_A,]
  data_B <- data[ix_B,]

  m_A <- glm(formula, data_A, family=binomial(link='logit'))
  yhat_A <- predict(m_A, data_A, type='link')
  g_A <- performance(yhat_A, m_A$y, ngroups = ngroups, ...)$gini

  m_B <- glm(formula, data_B, family=binomial(link='logit'))
  yhat_B <- predict(m_B, data_B, type='link')
  g_B <- performance(yhat_B, m_B$y, ngroups = ngroups, ...)$gini

  yhat_ALL <- c(yhat_A, yhat_B)
  y_ALL <- c(m_A$y, m_B$y)
  g_ALL <- performance(yhat_ALL, y_ALL, ngroups = ngroups, ...)$gini

  data.frame(
    variable = segvar,
    population = nrow(data),
    p_pob_A = nrow(data_A)/nrow(data),
    p_pob_B = nrow(data_B)/nrow(data),
    gini_TOT = g0,
    gini_A_B = g_ALL,
    gini_A = g_A,
    gini_B = g_B,
    p_pos_TOT = mean(m0$y),
    p_pos_A = mean(m_A$y),
    p_pos_B = mean(m_B$y),
    stringsAsFactors = F
  )
}

split_one(y ~ x + z, d, 'class1')


split.formula <- function(formula, data, segvars, ngroups = 100, ...){
  s <- lapply(segvars, function(v){
    split_one(formula, data, v, ngroups, ...)
  })
  tab <- do.call(rbind, s) %>%
    arrange(desc(gini_A_B)) %>%
    mutate(rank = row_number())
  out <- list(
    population = tab$population[1],
    p_pos_TOT = tab$p_pos_TOT[1],
    gini_TOT = tab$gini_TOT[1],
    table = tab[c('rank','variable','p_pob_A','p_pob_B','gini_A_B',
                  'gini_A','gini_B','p_pos_A','p_pos_B')]
  )
  class(out) <- 'segtree.split'
  out
}

split(y ~ x + z, d, c('class1','class2'))

# Useless for the moment
fork <- function(data, segvar){
  classes <- levels(as.factor(data[[segvar]]))[1:2]
  list(
    data_A = data[data[[segvar]] == classes[1]],
    data_B = data[data[[segvar]] == classes[1]]
  )
}

leaf <- function(segvars,
                 levels,
                 name=paste0('leaf_',as.numeric(Sys.time())),
                 splits=list()){
  out <- list(
    segvars = segvars,
    levels = levels,
    name = name,
    splits = splits
  )
  class(out) <- 'leaf'
  out
}

print.leaf <- function(ll, ...){
  if(is.null(ll$segvars)){
    def <- 'Root node'
  } else{
    def <- paste(paste(ll$segvars, ll$levels, sep=' = '),collapse=' > ')
  }
  cat(sprintf('----->>\nLeaf: %s\nDefinition:\n\t%s\nSplit:\n\n',
              ll$name, def))
  print(ll$splits, ...)
  cat('<<-----\n')
}
tt$leaves$root
tt2$leaves$L1_A
tt3$leaves$L1_A_L2_A

segtree <- function(formula, data, segvars, fast=FALSE, ...){
  if(fast){
    splits = list()
  } else{
    splits=split.formula(formula, data, segvars, ...)
  }
  out <- list(
    formula = formula,
    leaves = list('root'=leaf(segvars=NULL, levels=NULL, name='root',
                              splits=splits)),
    data = data,
    segvars = segvars,
    fast = fast
  )
  class(out) <- c('segtree')
  out
}

print.segtree <- function(tree, ...){
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ >>\nSegmentation Tree\n\nNumber of leaves: %d\nTarget: %s\nRegression Variables:\n\t%s\nAvailable segmentation variables:\n\t%s\nPopulation: %s\nLeaves:\n\n',
              length(tree$leaves),
              as.character(tree$formula[2]),
              as.character(tree$formula[3]),
              paste(tree$segvars, collapse = ', '),
              format(nrow(tree$d), scientific = F, big.mark = ',')))
  for(l in tree$leaves){
    print(l, global_pop = nrow(tree$data))
  }
  cat('<< ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~')
}

summary.segtree <- function(tree){
  out <- list(
    population = nrow(tree$data),
    leaves = length(tree$leaves),
    formula = tree$formula,
    segvars = tree$segvars
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
  cat(sprintf('~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ >>\nSegmentation Tree Summary\n\nNumber of leaves: %d\nTarget: %s\nRegression Variables:\n\t%s\nAvailable segmentation variables:\n\t%s\nPopulation: %s\nLeaf summary:\n\n',
              s$leaves,
              as.character(s$formula[2]),
              as.character(s$formula[3]),
              paste(s$segvars, collapse = ', '),
              format(s$population, scientific = F, big.mark = ',')))
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

  split.formula(tree$formula, dat, segvars)
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
    s$table[p_vars] <- apply(s$table[p_vars], 2, scales::percent)
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
    name = names[1]
  )
  leaf_B <- leaf(
    segvars = c(leaf$segvars, segvar),
    levels = c(leaf$levels, classes[2]),
    name = names[2]
  )
  new_leaves <- list(leaf_A,leaf_B)
  names(new_leaves) <- names[1:2]
  tree$leaves <- c(tree$leaves[-leaf_ix], new_leaves)
  if(!fast){
    tree$leaves[[names[1]]]$splits <- split.segtree(tree, names[1])
    tree$leaves[[names[2]]]$splits <- split.segtree(tree, names[2])
  }
  tree
}















