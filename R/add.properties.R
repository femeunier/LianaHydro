#' @name add.properties
#' @title add.properties
#' @author Félicien Meunier
#' @export
#' @description Add 2 list of models
#' @param models models
#' @param x some PLC

add.properties <- function(models, x = c(12,50,88)){

  Nmodels <- length(unlist(lapply(models,'[[','RMSE')))

  for (i in seq(1,Nmodels)){

    function.name <- names(models)[i]
    function.slope <- match.fun(paste0("slope.",function.name))
    function.param <- models[[i]][["best.params"]]
    function.invert <- match.fun(paste0("invert.",function.name))

    slopes <- do.call(function.slope,list(x,function.param[1],function.param[2]))
    invert <- do.call(function.invert,list(x,function.param[1],function.param[2]))

    names(x) <- paste0("PLC",x)
    names(slopes) <- paste0("ax",x)
    names(invert) <- paste0("P",x)

    models[[i]][["x"]] <- x
    models[[i]][["slopes"]] <- slopes
    models[[i]][["invert"]] <- invert

  }

  return(models)

}
