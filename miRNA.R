#' s3RNA
#' This stand-alone package allows for analysis of high-throughput, small RNA sequencing data of the pig (Sus scrofa) ovarian follicle used in Pan, Toms, et al. (2020)
#'
#' @export miRNA
#'
#' @return shiny application object
#'
#' @example \dontrun {mymiR <- miRNA()}
#'
#'


# return the microRNA dataset
miRNA <- function() {
 return  miRReads
}
