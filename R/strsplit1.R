#' Title
#'
#' @param x A char vector
#' @param split What to split on
#'
#' @return A char vector
#' @export
#'
#' @examples
#' x <- "alfa,bravo,charlie,delta"
#' strsplit1(x,split=",")
strsplit1 <- function(x, split) {
  strsplit(x, split = split)[[1]]
}
