
require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)

leaf <- function(segvars,
                 levels,
                 name=paste0('leaf_',as.numeric(Sys.time())),
                 splits=list(),
                 children=c(),
                 terminal=FALSE){
  out <- list(
    segvars = segvars,
    levels = levels,
    name = name,
    splits = splits,
    children = children,
    terminal = terminal
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
  terminal <- ifelse(ll$terminal, 'Yes', 'No')
  cat(sprintf('----->>\nLeaf: %s\nTerminal: %s\nDefinition:\n\t%s\nChildren: %s\nSplit:\n\n',
              ll$name, terminal, def, paste(ll$children, collapse=', ')))
  print(ll$splits, ...)
  cat('<<-----\n')
}
