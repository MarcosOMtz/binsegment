
require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)

split_one_ <- function(formula, data, segvar, ngroups = 100, ...){
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
    p_pob_A = nrow(data_A)/nrow(data),
    p_pob_B = nrow(data_B)/nrow(data),
    gini_A_B = g_ALL,
    gini_A = g_A,
    gini_B = g_B,
    p_pos_A = mean(m_A$y),
    p_pos_B = mean(m_B$y),
    stringsAsFactors = F
  )
}

split_.formula <- function(formula, data, segvars, ngroups = 100, ...){
  m0 <- glm(formula, data, family=binomial(link='logit'))
  yhat0 <- predict(m0, data, type='link')
  g0 <- performance(yhat0, m0$y, ngroups = ngroups, ...)$gini

  s <- lapply(segvars, function(v){
    split_one_(formula, data, v, ngroups, ...)
  })
  tab <- do.call(rbind, s) %>%
    arrange(desc(gini_A_B)) %>%
    mutate(rank = row_number())
  out <- list(
    population = nrow(data),
    p_pos_TOT = mean(as.numeric(as.character(data[[as.character(formula)[2]]]))),
    gini_TOT = g0,
    table = tab[c('rank','variable','p_pob_A','p_pob_B','gini_A_B',
                  'gini_A','gini_B','p_pos_A','p_pos_B')]
  )
  class(out) <- 'segtree.split'
  out
}
