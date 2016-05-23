
require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)

percent <- function(x, digits=0){
  fmt <- sprintf('%%.%df', digits)
  sprintf(paste0(fmt,'%%'), 100*x)
}

split_one_ <- function(formula, data, segvar, ngroups = 100, ...){
  classes <- levels(as.factor(data[[segvar]]))[1:2]
  ix_A <- data[[segvar]] == classes[1]
  ix_B <- data[[segvar]] == classes[2]
  if(length(unique(data[[segvar]])) < 2 || sum(ix_A) == 0 || sum(ix_B) == 0){
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

  A <- classes[1]
  B <- classes[2]
  short_A <- ifelse(nchar(A) > 10, paste0(substr(A, 1, 7), '...'), A)
  short_B <- ifelse(nchar(B) > 10, paste0(substr(B, 1, 7), '...'), B)
  data.frame(
    variable = segvar,
    A = short_A,
    B = short_B,
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
  if(all(sapply(s, is.null))){
    tab <- data.frame(rank=0,variable='N/A',A='N/A',B='N/A',p_pob_A=0,p_pob_B=0,gini_A_B=0,gini_A=0,gini_B=0,p_pos_A=0,p_pos_B=0)
  } else{
    tab <- do.call(rbind, s) %>%
      arrange(desc(gini_A_B)) %>%
      mutate(rank = row_number())
    tab <- tab[c('rank','variable','A','B','p_pob_A','p_pob_B','gini_A_B',
                 'gini_A','gini_B','p_pos_A','p_pos_B')]
  }

  out <- list(
    population = nrow(data),
    p_pos_TOT = mean(as.numeric(as.character(data[[as.character(formula)[2]]]))),
    gini_TOT = g0,
    table = tab
  )
  class(out) <- 'segtree.split'
  out
}

leaf_index <- function(leaf, data){
  aux <- sapply(data[leaf$segvars],
                function(x) trimws(as.character(x), which = "both"))
  apply(aux, 1, function(x){
    all(x == leaf$levels)
  })
}

split_data <- function(tree, index.only = T, newdata=NULL, ...){
  if(is.null(newdata)){
    newdata <- tree$data
  }
  newdata <- as.data.frame(newdata)
  leaves <- tree$leaves[sapply(tree$leaves, function(x) x$terminal)]
  codes <- as.data.frame(lapply(1:length(leaves), function(i){
    i*leaf_index(leaves[[i]], newdata)
  })) %>%
    apply(1, sum)
  if(any(codes == 0)){
    nms <- c('N/A', names(leaves))
    codes <- codes + 1
  } else{
    nms <- names(leaves)
  }
  codes <- nms[codes]

  if(index.only){
    return(codes)
  } else{
    return(split.data.frame(newdata, as.factor(codes), ...))
  }
}

structure.segtree <- function(tree, leaf_name = 'root'){
  leaf <- tree$leaves[[leaf_name]]
  depth <- length(leaf$segvars)
  if(leaf$terminal){
    if(depth == 0){
      def <- sprintf('(%s)', leaf$name)
    } else{
      def <- sprintf('> %s = %s (%s)',
                   leaf$segvars[depth], leaf$levels[depth], leaf$name)
    }
    s <- list(name=leaf_name,
              def=def,
              details=list(
                p_population=leaf$splits$population/nrow(tree$data),
                p_pos=leaf$splits$p_pos_TOT,
                gini=leaf$splits$gini_TOT
              ),
              children=NULL)
  } else{
    if(depth == 0){
      def <- sprintf('(%s)', leaf$name)
    } else{
      def <- sprintf('> %s = %s (%s)', leaf$segvars[depth], leaf$levels[depth], leaf$name)
    }
    s <- list(name=leaf_name,
              def=def,
              details=list(
                p_population=leaf$splits$population/nrow(tree$data),
                p_pos=leaf$splits$p_pos_TOT,
                gini=leaf$splits$gini_TOT
              ),
              children=lapply(leaf$children, function(l) structure.segtree(tree, l)))
  }
  class(s) <- 'structure.segtree'
  s
}

print.structure.segtree <- function(s, details = TRUE, prefix = '', level = 0){
  det <- ''
  if(details){
    if(level <= 0){
      cat('[% pop, p(Y = 1), Gini]\n\n')
    }
    if(is.null(s$children)){
      det <- sprintf(' --> [%s, %s, %.1f]',
                     percent(s$details$p_population, 1),
                     percent(s$details$p_pos, 1),
                     100*s$details$gini)
    }
  } else{
    if(level <= 0){
      cat('\n')
    }
  }
  cat(sprintf('%s%s %s\n', prefix, s$def, det))
  if(!is.null(s$children)){
    for(c in s$children){
      print(c, details=details, prefix=paste(prefix, '  '), level=level+1)
    }
  }
}


