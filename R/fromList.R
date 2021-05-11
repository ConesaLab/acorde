# MODIFIED fromList() FROM UpSetR() PACKAGE THAT KEEPS ROWNAMES

#' @export
fromList <- function (input) {
  # original fromList()
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # preserve rownames
  row.names(data) <- elements
  return(data)
}
