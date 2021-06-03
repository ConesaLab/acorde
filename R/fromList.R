#' @title Obtaining intersections from a list of named vectors
#'
#' @description Modified version of \code{\link[UpSetR]{fromList}}
#' (\code{UpSetR} package). Transforms a list of named vectors to a data frame
#' compatible with the \code{\link[UpSetR]{upset}} function.
#'
#' @param list A list of named vectors, each including one group of elements of
#' interest.
#'
#' @return A data.frame with the (unique) elements within the list as row names
#' and one column per group. Values of 0 or 1 will be shown if the element is
#' present or absent from the group, respectively.
#'
#' @references
#'
#' The \code{fromList()} function is originally included in the
#' \code{\link[UpSetR:fromList]{UpSetR}} package:
#'
#' \insertRef{Conway2017}{acorde}
#'
#' The present modifications to the fromList() function were obtained from
#' user @@docmanny's contribution in this
#' \href{https://github.com/hms-dbmi/UpSetR/issues/85#issuecomment-327900647}{GitHub issue}.
#'
#' @export
fromList <- function (list) {
  # original fromList()
  elements <- unique(unlist(list))
  data <- unlist(lapply(list, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(list), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(list)
  # preserve rownames
  row.names(data) <- elements
  return(data)
}
