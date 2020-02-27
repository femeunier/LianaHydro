#' @name default.model.params
#' @title default.model.params
#' @author Félicien Meunier
#' @export
#' @description Returns default model param
#' @param fun.name Function name

default.model.params <- function(fun.name){

  switch(fun.name,
         weibull = {
           list(2,2)
         },
         sigmoidal = {
           list(2,-2)
         },
         polynomial2 = {
           list(2,-2)
         },
         polynomial = {
           list(2,-2)
         },
         stop("Enter something that switches me!"))
  }
