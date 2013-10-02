\name{var2rV}
\alias{var2rV}
\title{Variance component correlations}
\usage{
  var2rV(V)
}
\description{
  Turn an array of variances/covariance matrices into
  correlations
}
\examples{
data(BTdata)
m <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
V <- reshapeVCV(m$VCV)
rV <- var2rV(V)
}

