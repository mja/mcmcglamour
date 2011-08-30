
match.length <- function(x) attr(x, 'match.length')

match.random <- function(x) {

        r <- regexpr("(at.level\\(\\S+, \\d+\\)\\.)??\\w+$", x, perl=TRUE) # match from the last : to the end of the line
        random <- substr(x, r[1], r[1] + match.length(r)[1] - 1)

        return(random)
 }

 match.response <- function(x) {

         r <- regexpr("^[a-z]+", x, perl=TRUE)
         response <- substr(x, r[1], r[1] + match.length(r)[1] - 1)

         return(response)
 }

# Turn a two dimensional nitt × (v × p × p) matrix of (co)variances
# into a nitt × p × p × v dimensional array
reshape.VCV <- function(VCV) {

  # extract random effect and response terms from the VCV dimension names
  random <- unique(sapply(dimnames(VCV)[[2]], match.random))
  response <- unique(sapply(dimnames(VCV)[[2]], match.response))

  nitt <- dim(VCV)[1] # number of iterations
  V <- array(as.vector(VCV), dim=c(nitt, length(response), length(response), length(random)), dimnames=list(itt=as.character(1:nitt), trait=response, trait=response, V=random))

  return(V)
}

# Proportional variance compnents, v^2
var2v2 <- function(V) {

  require(plyr)

  vcv.diag <- plyr::aaply(V, .margins=c(1,4), .fun=diag)

  v2 <- plyr::aaply(vcv.diag, .margins=c(1,3), .fun=function(V) V / sum(V))

  return(v2)
}

# Turn into correlations
var2rV <- function(V) {

   require(plyr)

   rV <- plyr::aaply(V, .margins=c(1,4), .fun=cov2cor)

   return(rV)
}
