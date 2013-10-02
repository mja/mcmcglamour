\name{reshapeVCV}
\alias{reshapeVCV}
\title{Reshape MCMCglmm VCV matrix into an array}
\usage{
  reshapeVCV(VCV)
}
\description{
  Turn a two dimensional nitt × (v × p × p) matrix of
  (co)variances into a nitt × p × p × v dimensional array
}
\examples{
data(BTdata)
m <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
V <- reshapeVCV(m$VCV)
}

