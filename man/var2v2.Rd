\name{var2v2}
\alias{var2v2}
\title{Proportional variance components}
\usage{
  var2v2(V)
}
\description{
  Proportional variance components, v^2, from an array of
  variances and covariances
}
\examples{
data(BTdata)
m <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
V <- reshapeVCV(m$VCV)
v2 <- var2v2(V)
}

