#' @name extract.model
#' @title extract.model
#' @author Félicien Meunier
#' @export
#' @description Returns model
#' @param all.models all.models
#' @param pos pos


extract.model <- function(all.models,
                          pos){

  return(all.models[[pos]])
}
