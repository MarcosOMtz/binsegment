

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
         class2 = factor(ifelse(runif(ng+nb) < abs(x*z/max(x*z)), 'C', 'D')))

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
  wr0 <- wroc(yhat0, m0$y, ngroups = 100, ...)

  classes <- levels(data[[segvar]])[1:2]

  data_A <- data[data[[segvar]] == classes[1],]
  m_A <- glm(formula, data_A, family=binomial(link='logit'))
  yhat_A <- predict(m_A, data_A, type='link')
  wr_A <- wroc(yhat_A, m_A$y, ngroups = 100, ...)

  data_B <- data[data[[segvar]] == classes[2],]
  m_B <- glm(formula, data_B, family=binomial(link='logit'))
  yhat_B <- predict(m_B, data_B, type='link')
  wr_B <- wroc(yhat_B, m_B$y, ngroups = 100, ...)

  yhat_ALL <- c(yhat_A, yhat_B)
  y_ALL <- c(m_A$y, m_B$y)
  wr_ALL <- wroc(yhat_ALL, y_ALL, ngroups = 100, ...)

  data.frame(
    variable = segvar,
    poblacion = nrow(data),
    p_pob_A = nrow(data_A)/nrow(data),
    p_pob_B = nrow(data_B)/nrow(data),
    gini_TOT = performance(wr0)$gini,
    gini_A_B = performance(wr_ALL)$gini,
    gini_A = performance(wr_A)$gini,
    gini_B = performance(wr_B)$gini,
    stringsAsFactors = F
  )
}

split_one(y ~ x + z, d, 'class1')


split <- function(formula, data, segvars, ngroups = 100, ...){
  s <- lapply(segvars, function(v){
    split_one(formula, data, v, ngroups, ...)
  })
  do.call(rbind, s)
}

split(y ~ x + z, d, c('class1','class2'))












