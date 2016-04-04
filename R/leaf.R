
require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)

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
