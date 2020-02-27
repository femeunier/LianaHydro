#' @name fill.NA.data
#' @title fill.NA.data
#' @author Félicien Meunier
#' @export
#' @description Returns full data.frame
#' @param data.in data.in

fill.NA.data <- function(data.in = data){

  if(!all(!is.na(data.in[1,]))) stop("First line not full of values")

  data.out <- data.in
  cols <- colnames(data.in)
  Nrow <- nrow(data.in)

  for (col in cols){
    for (irow in seq(1,Nrow)){
      if(is.na(data.in[irow,col])){
        data.out[irow,col] <-  data.out[irow-1,col]
      }
    }
  }

  return(data.out)
}
