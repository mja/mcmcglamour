#' Read in a vector of RDS as a list
#'
#' @export
#' @examples
#' read_rdss(c('one.rds', 'two.rds'))
read_rdss <- function(filepaths) plyr::alply(filepaths, 1, readRDS)
