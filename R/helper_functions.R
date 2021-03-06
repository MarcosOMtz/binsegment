
percent <- function(x, digits=0){
  fmt <- sprintf('%%.%df', digits)
  sprintf(paste0(fmt,'%%'), 100*x)
}

split_one_ <- function(formula, data, segvar, ngroups = 100, ...){
  fac <- as.factor(data[[segvar]])
  unq <- unique(fac)
  unq <- unq[!is.na(unq)]
  if(length(unq) < 2){
    warning(sprintf('Segmentation variable %s has less than two unique values available. Skipping...', segvar))
    return(NULL)
  } else if(length(unq) > 2){
    warning(sprintf('Segmentation variable %s has more than two unique values. Using only the first two: "%s", "%s"', segvar, as.character(levels(fac)[1]), as.character(levels(fac)[2])))
  }
  classes <- levels(fac)[1:2]

  ix_A <- data[[segvar]] == classes[1]
  ix_B <- data[[segvar]] == classes[2]
  if(any(is.na(ix_A)) | any(is.na(ix_B))){
    warning(sprintf('Segmentation variable %s NA values. Removing...', segvar))
    ix_A <- ix_A[!is.na(ix_A)]
    ix_B <- ix_B[!is.na(ix_B)]
  }
  data_A <- data[ix_A,]
  data_B <- data[ix_B,]

  m_A <- glm(formula, data_A, family=binomial(link='logit'))
  yhat_A <- predict(m_A, data_A, type='link')
  g_A <- optimbucket::performance(yhat_A, m_A$y, ngroups = ngroups, ...)$gini

  m_B <- glm(formula, data_B, family=binomial(link='logit'))
  yhat_B <- predict(m_B, data_B, type='link')
  g_B <- optimbucket::performance(yhat_B, m_B$y, ngroups = ngroups, ...)$gini

  yhat_ALL <- c(yhat_A, yhat_B)
  y_ALL <- c(m_A$y, m_B$y)
  g_ALL <- optimbucket::performance(yhat_ALL, y_ALL, ngroups = ngroups, ...)$gini

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
  g0 <- optimbucket::performance(yhat0, m0$y, ngroups = ngroups, ...)$gini

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

#' @rdname split_data
#' @export
leaf_index <- function(leaf, data){
  aux <- sapply(data[leaf$segvars],
                function(x) trimws(as.character(x), which = "both"))
  apply(aux, 1, function(x){
    all(x == leaf$levels)
  })
}

#' Apply Segmentation to Data
#'
#' Returns a factor with the terminal node each observation corresponds to.
#' Alternatively, returns a list of \code{data.frame}s corresponding to each
#' terminal node. This function is similar to some predict methods.
#'
#' @param tree tree An object of class \code{segtree}
#' @param index.only Whether to return the segmentation variable (if
#'   \code{TRUE}, the default) or to actually split the data
#' @param newdata A \code{data.frame} to be split. If \code{NULL}, then the
#'   tree's data is used
#' @param leaf An object of class \code{leaf}
#' @param data [\code{leaf_index} only] A \code{data.frame} to be segmented
#' @details The function \code{split_data} a convenient wrapper of
#'   \code{leaf_index}. \code{leaf_index} simply returns a boolean vector of
#'   whether each observation is in the given leaf or not.
#' @export
split_data <- function(tree, index.only = T, newdata=NULL){
  if(is.null(newdata)){
    message("No data provided. Splitting the tree's data.")
    newdata <- tree$data
  } else if(!is.matrix(newdata) && !is.data.frame(newdata)){
    stop('newdata must be a matrix or data.frame, not a vector.')
  }
  newdata <- as.data.frame(newdata)
  actual_segvars <- unique(unlist(lapply(tree$leaves, function(x) x$segvars)))
  if(!(all(actual_segvars %in% names(newdata)))){
    mis <- setdiff(actual_segvars, names(newdata))
    stop(sprintf('The following variables are not present in newdata: %s',
         paste(mis, collapse = ', ')))
  }
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
    return(split.data.frame(x=newdata, f=as.factor(codes)))
  }
}


#' Structure of a Segmentation Tree
#'
#' Shows a tree-like data structure with the dependences between the nodes of a
#' \code{segtree}.
#' @param tree tree An object of class \code{segtree}
#' @param leaf_name The name of a leaf in the tree
#' @param details The level of details to show in the tree structure. 0 means no
#'   details, 1 means only on terminal nodes, 2 means on all nodes and >= 3
#'   means on all nodes plus extra information
#' @return An object of class \code{structure.segtree}, which draws the tree's
#'   structure in a text graph
#' @export
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
                population = leaf$splits$population,
                tot_population = nrow(tree$data),
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
                population = leaf$splits$population,
                tot_population = nrow(tree$data),
                p_population=leaf$splits$population/nrow(tree$data),
                p_pos=leaf$splits$p_pos_TOT,
                gini=leaf$splits$gini_TOT
              ),
              children=lapply(leaf$children, function(l) structure.segtree(tree, l)))
  }
  class(s) <- 'structure.segtree'
  s
}

#' @rdname structure.segtree
#' @export
print.structure.segtree <- function(s, details = c(0, 1, 2, 3), prefix = '', level = 0){
  det <- ''
  if(level <= 0){
    if(details[1] >= 1){
      if(details[1] <= 2){
        cat('[% pop, p(Y = 1), Gini]\n\n')
      } else if(details[1] >=3){
        cat('[pop (% pop), bads (p(Y = 1)), Gini]\n\n')
      }
    }
  }
  if(details[1] >= 1){
    if(details[1] <=2){
      if(details[1] >= 2 | is.null(s$children)){
        det <- sprintf(' --> [%s, %s, %.1f]',
                       percent(s$details$p_population, 1),
                       percent(s$details$p_pos, 1),
                       100*s$details$gini)
      }
    } else if(details[1] >= 3){
        det <- sprintf(' --> [%s (%s), %s (%s), %.1f]',
                       format(s$details$population, big.mark = ','),
                       percent(s$details$p_population, 1),
                       format(round(s$details$p_pos*s$details$population), big.mark = ','),
                       percent(s$details$p_pos, 1),
                       100*s$details$gini)
    }

  }
  cat(sprintf('%s%s %s\n', prefix, s$def, det))
  if(!is.null(s$children)){
    for(c in s$children){
      print(c, details=details, prefix=paste(prefix, '  '), level=level+1)
    }
  }
}


