

require(dplyr)
require(tidyr)
require(ggplot2)
library(optimbucket)

#### Data 1
ng <- 10000
nb <- 1000

goods <- data.frame(
  x = round(rnorm(ng, mean = -2),2),
  z = round(runif(ng, -3, 3)),
  y = 0
)
bads <- data.frame(
  x = rnorm(nb, 0, 1),
  z = round(runif(nb, 0, 1)),
  y = 1
)

d <- rbind(goods, bads) %>%
  mutate(y = factor(y),
         class1 = factor(ifelse(runif(ng+nb) < abs(x/max(x)), 'A', 'B')),
         class2 = factor(ifelse(runif(ng+nb) < abs(x*z/max(x*z)), 'C', 'D')),
         class3 = factor(ifelse(runif(ng+nb) < abs((x+z)/max(x+z)), 'X', 'Z')))

ggplot(d, aes(x, fill=y)) +
  geom_density(alpha=0.5) +
  facet_wrap(~class1)

ggplot(d, aes(z, fill=y)) +
  geom_density(alpha=0.5) +
  facet_wrap(~class1)


m <- glm(y ~ x + z, data = d, family = 'binomial')
s <- summary(m)

d2 <- cbind(d, prob = m$fitted.values)
ggplot(d2, aes(z, prob, fill=y)) +
  geom_point() +
  facet_wrap(~class1)


glms <- function(formula, data, segvar, family=binomial(link='logit'), ...){
  if(!is.factor(data[[segvar]])){
    warning('Segmentation variable must be a factor. Coercing to factor')
    data[[segvar]] <- as.factor(data[[segvar]])
  }
  classes <- levels(data[[segvar]])[1:2]
  models <- list()
  models[[1]] <- glm(formula,
                     data=data[data[[segvar]] == classes[1],],
                     family=family,
                     ...)
  models[[2]] <- glm(formula,
                     data=data[data[[segvar]] == classes[2],],
                     family=family,
                     ...)
  covars <- list()
  covars[[1]] <- vcov(models[[1]])
  covars[[2]] <- vcov(models[[2]])

  names(models) <- names(covars) <- classes

  out <- list(
    classes = classes,
    models = models,
    covars = covars,
    formula = formula,
    call = match.call()
  )
  class(out) <- c('list', 'glms')
  out
}

m <- glms(y ~ x + z, data = d, 'class1')


logistic <- function(z){
  as.vector(1/(1 + exp(-z)))
}

logistic_derivative <- function(z){
  s <- logistic(z)
  s*(1-s)
}

plot(logistic, xlim=c(-4,4))
plot(logistic_derivative, xlim=c(-4,4))

var_probs <- function(x, beta, var_beta, add_ones = (length(beta) = ncol(x) + 1)){
  if(is.vector(x)){
    x <- matrix(x, nrow = 1)
  } else if(is.data.frame(x)){
    x <- as.matrix(x)
  }
  if(add_ones){
    x <- cbind(1, x)
  }
  xf <- logistic_derivative(x %*% beta)*x
  apply(xf, 1, function(v){
    v %*% var_beta %*% v
  })
  # xf %*% var_beta %*% t(xf)
}


var_probs(d[1:10,c('x','z')], m$models$A$coefficients, m$covars$A)


segment_test <- function(models, data, significance = 0.05){
  x <- model.matrix(models$formula, data)
  probs <- list()
  probs[[1]] <- logistic(x %*% models$models[[1]]$coefficients)
  probs[[2]] <- logistic(x %*% models$models[[2]]$coefficients)
  vars <- list()
  vars[[1]] <- var_probs(x,
                         models$models[[1]]$coefficients,
                         models$covars[[1]],
                         add_ones = F)
  vars[[2]] <- var_probs(x,
                         models$models[[2]]$coefficients,
                         models$covars[[2]],
                         add_ones = F)
  z <- (probs[[1]] - probs[[2]])/sqrt(vars[[1]] + vars[[2]])
  pval <- pnorm(z, 0, 1)
  pval <- 2*pmin(pval, 1-pval)
  reject_null <- (pval <= significance)

  out <- list(
    statistic = z,
    p_value = pval,
    reject_null = reject_null,
    num_rejections = sum(reject_null),
    prob_reject = mean(reject_null),
    mean_p_value = mean(pval),
    sd_p_value = sd(pval)
  )
  class(out) <- c('list', 'segment_test')
  out
}

a <- segment_test(m, d, 0.05)


split <- function(formula, data, segvars, significance = 0.05,
                  family = binomial(link='logit'), ...){
  out <- lapply(segvars, function(v){
    gs <- glms(formula, data, v, family, ...)
    test <- segment_test(gs, data, significance)
    temp <- as.data.frame(
      test[c('num_rejections', 'prob_reject', 'mean_p_value', 'sd_p_value')]
    )
    cbind(variable=v, temp)
  })
  do.call(rbind, out)
}

split(y~x+z, d, c('class1','class2'))





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
    poblacion = nrow(data),
    p_pob_A = nrow(data_A)/nrow(data),
    p_pob_B = nrow(data_B)/nrow(data),
    gini_TOT = g0,
    gini_A_B = g_ALL,
    gini_A = g_A,
    gini_B = g_B,
    tm_TOT = mean(m0$y),
    tm_A = mean(m_A$y),
    tm_B = mean(m_B$y),
    stringsAsFactors = F
  )
}

split_one(y ~ x + z, d, 'class1')


split.formula <- function(formula, data, segvars, ngroups = 100, ...){
  s <- lapply(segvars, function(v){
    split_one(formula, data, v, ngroups, ...)
  })
  do.call(rbind, s)
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
    segvars = segvars
  )
  class(out) <- c('segtree')
  out
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

fork.segtree <- function(tree, leaf, segvar, names, fast=FALSE){
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

tt <- segtree(y ~ x + z, (d), c('class1','class2','class3'), fast=F)
#split.segtree(tt, 'root')
tt2 <- fork.segtree(tt, 'root', 'class3', c('L1_A', 'L1_B'), fast = F)
split.segtree(tt2, 'L1_A')
tt3 <- fork.segtree(tt2, 'L1_A', 'class2', c('L1_A_L2_A', 'L1_A_L2_B'))















