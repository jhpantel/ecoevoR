#' round a number up to the nearest power of 10
#'
#' @param x the number to round
#'
#' @return the rounded value
#' @export
#'
#' @examples
#' round_up(9)
#' round_up(1020)
round_up <- function(x){
  10^ceiling(log10(x))
}
