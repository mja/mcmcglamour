\name{merge_chains}
\alias{merge_chains}
\title{Merge a chained-linked MCMCglmm model}
\usage{
  merge_chains(model)
}
\description{
  Turns a linked MCMCglmm object produced by
  \code{\link{link_chains}} into an object with
  concatenated \code{Sol} and \code{VCV} matrices.
}
\examples{
data(BTdata)
m1 <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
m2 <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
m <- link_chains(list(m1, m2))
bt_model <- merge_chains(m)
summary(bt_model)
}

