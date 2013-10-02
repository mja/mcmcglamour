\name{link_chains}
\alias{link_chains}
\title{Link MCMCglmm objects as a list}
\usage{
  link_chains(chains)
}
\description{
  Links together multiple runs (chains) of an MCMCglmm
  model. It turns the Sol and VCV slots into mcmc.list
  objects
}
\examples{
data(BTdata)
m1 <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
m2 <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
m <- link_chains(list(m1, m2))
plot(m$VCV)
}

