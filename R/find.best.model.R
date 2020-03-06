#' @name find.best.model
#' @title find.best.model
#' @author Félicien Meunier
#' @export
#' @description Returns best model
#' @param data.test data
#' @param function.names function.names

find.best.model <- function(All.models){

  best.modelpos <- which.min(unlist(lapply(All.models,'[[','RMSE')))
  best.model <- extract.model(All.models,best.modelpos)

  return(list(best.model))
}
