#' @name plot.models
#' @title plot.models
#' @author Félicien Meunier
#' @export
#' @description plot models
#' @param models models

plot.models <- function(models,add =TRUE,col = 'black',highlight = FALSE){

  Nmodels <- length(unlist(lapply(models,'[[','RMSE')))
  best.modelpos <- which.min(unlist(lapply(models,'[[','RMSE')))

  if (length(col)<Nmodels) col<- rep(col,Nmodels)

  if (highlight) col[best.modelpos] <- 'red'

  for (i in seq(1,Nmodels)){
    if (!add & i == 1){
      plot(models[[i]][["psi.all"]],models[[i]][["PLC.predict.all"]],col = col[i])

    } else{
      lines(models[[i]][["psi.all"]],models[[i]][["PLC.predict.all"]],col= col[i])
    }
  }


}
