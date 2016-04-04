require(dplyr)
require(tidyr)
require(ggplot2)
require(optimbucket)

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


# Usage

tt <- segtree(y ~ x + z, (d), c('class1','class2','class3'), fast=F)
# split.segtree(tt, 'root')
tt$leaves$root$splits
tt2 <- fork.segtree(tt, 'root', 'class3', c('L1_A', 'L1_B'), fast = F)
# split.segtree(tt2, 'L1_A')
tt3 <- fork.segtree(tt2, 'L1_A', 'class2', c('L1_A_L2_A', 'L1_A_L2_B'))

############
# Function tests
split_one(y ~ x + z, d, 'class1')
split(y ~ x + z, d, c('class1','class2'))
