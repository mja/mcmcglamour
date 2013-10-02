
match_length <- function(x) attr(x, 'match.length')

match_random <- function(x) {

        r <- regexpr("(at.level\\(\\S+, \\d+\\)\\.)??\\w+$", x, perl=TRUE) # match from the last : to the end of the line
        random <- substr(x, r[1], r[1] + match_length(r)[1] - 1)

        return(random)
 }

 match_response <- function(x) {

         r <- regexpr("^[a-zA-z_0-9]+", x, perl=TRUE)
         response <- substr(x, r[1], r[1] + match_length(r)[1] - 1)

         return(response)
 }

#' Reshape MCMCglmm VCV matrix into an array
#'
#' Turn a two dimensional nitt × (v × p × p) matrix of (co)variances
#' into a nitt × p × p × v dimensional array
#'
#' @export
#' @examples
#' data(BTdata)
#' m <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
#' V <- reshapeVCV(m$VCV)
reshapeVCV <- function(VCV) {

  # extract random effect and response terms from the VCV dimension names
  random <- unique(sapply(dimnames(VCV)[[2]], match_random))
  response <- unique(sapply(dimnames(VCV)[[2]], match_response))

  nitt <- dim(VCV)[1] # number of iterations
  V <- array(as.vector(VCV), dim=c(nitt, length(response), length(response), length(random)), dimnames=list(itt=as.character(1:nitt), trait1=response, trait2=response, V=random))

  return(V)
}

reshape.VCV <- reshapeVCV

#' Proportional variance components
#'
#' Proportional variance components, v^2, from an array of variances
#' and covariances
#' @export
#' @examples
#' data(BTdata)
#' m <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
#' V <- reshapeVCV(m$VCV)
#' v2 <- var2v2(V)
var2v2 <- function(V) {

  require(plyr)

  vcv_diag <- plyr::aaply(V, .margins=c(1,4), .fun=diag)

  v2 <- plyr::aaply(vcv_diag, .margins=c(1,3), .fun=function(V) V / sum(V))

  return(v2)
}

#' Variance component correlations
#'
#' Turn an array of variances/covariance matrices into correlations
#' @export
#' @examples
#' data(BTdata)
#' m <- MCMCglmm(cbind(tarsus, back) ~ trait + trait:sex + trait:hatchdate - 1, random=~us(trait):dam, rcov=~us(trait):units, data=BTdata, family=c('gaussian', 'gaussian'))
#' V <- reshapeVCV(m$VCV)
#' rV <- var2rV(V)
# Turn into correlations
var2rV <- function(V) {

   require(plyr)

   rV <- plyr::aaply(V, .margins=c(1,4), .fun=cov2cor)

   return(rV)
}
