
require(testthat)
require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)
require(binsegment)

set.seed(1234)

data_example <- function(ng = 10000, nb = 1000){
  n <- ng + nb
  goods <- data.frame(
    x = round(rnorm(ng, mean = -2),2),
    z = round(runif(ng, -3, 3)),
    w = 3,
    a = NA,
    b = NaN,
    c = Inf,
    d = ifelse(runif(ng) < 0.1, NA, runif(ng)),
    y = 0
  )
  bads <- data.frame(
    x = rnorm(nb, 0, 1),
    z = round(runif(nb, 0, 1)),
    w = 3,
    a = NA,
    b = NaN,
    c = Inf,
    d = ifelse(runif(nb) < 0.1, NA, runif(nb)),
    y = 1
  )

  d <- rbind(goods, bads) %>%
    mutate(y = factor(y),
           class1 = sample(c('A','B'), size=n, replace = T),
           class2 = sample(1:2, size=n, replace = T),
           class3 = sample(c(TRUE, FALSE), size=n, replace = T),
           class4 = factor(class1),
           class5 = 'U',
           class6 = sample(LETTERS, size=n, replace = T),
           class7 = ifelse(runif(n) < 0.5, class6, NA),
           class8 = ifelse(runif(n) < 0.5, class2, NA))
  d
}

data_example(10,5)


### Useful everywhere
segvars <- paste('class', 1:8, sep='')

### segtree

test_that('segtree can handle invalid inputs', {
  d <- data_example(1000, 300)
  expect_error(segtree('y ~ x', data = d, segvars = segvars[1:4])) #
  expect_error(segtree(y ~ x, data = d[[1]], segvars = segvars[1:4]))
  # expect_error(segtree(y ~ x, data = d, segvars = c('x','a'))) #

  expect_error(segtree(y ~ x, data = d, segvars = c('hola','z')))
  expect_error(segtree(hola ~ x, data = d, segvars = c('x','a')))
})

test_that('segtree can handle bad regression variables', {
  d <- data_example(1000, 300)
  test_regvar <- function(v){
    f <- as.formula(paste('y', v, sep='~'))
    segtree(f, data=d, segvars = segvars[1:4])
  }

  expect_is(test_regvar('x'), 'segtree')
  expect_is(test_regvar('z'), 'segtree')
  expect_is(test_regvar('w'), 'segtree')
  expect_error(test_regvar('a'))
  expect_error(test_regvar('b'))
  expect_error(test_regvar('c'))
  expect_error(test_regvar('d'))
})

test_that('segtree can handle bad segmentation variables', {
  d <- data_example(1000, 300)
  test_segvar <- function(v){
    segtree(y ~ x, data=d, segvars = v)
  }

  expect_is(test_segvar('class1'), 'segtree')
  expect_is(test_segvar('class2'), 'segtree')
  expect_is(test_segvar('class3'), 'segtree')
  expect_is(test_segvar('class4'), 'segtree')
  expect_warning(test_segvar('class5'), 'one level|less than two')
  expect_warning(test_segvar('class6'), 'more than two')
  expect_warning(test_segvar('class7'), 'NA')
  expect_warning(test_segvar('class8'), 'NA')
})

test_that('fork can handle a variety of inputs', {
  d <- data_example(1000, 300)
  tr <- segtree(y ~ x + z, data=d, segvars = segvars[1:4])

  expect_error(fork(1:10, 'root', segvars[1], c('A','B')))
  expect_error(fork(tr, 'i_dont_exist_node', segvars[1], c('A','B')))
  expect_error(fork(tr, 'root', 'i_dont_exist_segvar', c('A','B')))
  expect_error(fork(tr, 'root', 'z', c('A','B'))) # isnt in segvars and too many levels
})

test_that('fork can set new node names automatically', {
  d <- data_example(1000, 300)
  tr <- segtree(y ~ x + z, data=d, segvars = segvars[1:4])

  expect_warning(fork(tr, 'root', segvars[1], 1), 'automatic|names') # a name is missing
  expect_is(fork(tr, 'root', segvars[1], c('a','b','c')), 'segtree')
  expect_message(fork(tr, 'root', segvars[1]))
})





















