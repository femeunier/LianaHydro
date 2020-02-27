#' @name opt.data
#' @title opt.data
#' @author Félicien Meunier
#' @export
#' @description Returns list
#' @param data data
#' @param function.names function.names


opt.data <-
  function(data = data,
         function.names = c("weibull","sigmoidal","polynomial","polynomial2")){

  function2test <- lapply(as.list(function.names),match.fun)

  models <- list()

  psi <- data$psi
  psi_all <- seq(min(psi),max(psi),length.out = 10000)

  for (i in seq(1,length(function2test))){

    models[[function.names[i]]] <- list()
    default <- default.model.params(function.names[i])
    m <- nlsLM(data = data,
               PLC ~ do.call(function2test[[i]],list(psi,param1,param2)),
               start=list(param1 = default[[1]], param2 = default[[2]]), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                                               printEval = TRUE, warnOnly = TRUE))
    models[[function.names[i]]][["m"]] <- m

    models[[function.names[i]]][["psi"]] <- psi
    models[[function.names[i]]][["psi.all"]] <- psi_all

    best.params <- coef(m)
    models[[function.names[i]]][["best.params"]] <- best.params

    RMSE <- do.call(match.fun(paste0(function.names[i],".comp")),
                    list(data,
                         best.params[1],
                         best.params[2]))$RMSE

    models[[function.names[i]]][["RMSE"]]  <- RMSE

    models[[function.names[i]]][["PLC.predict"]] <-
      do.call(match.fun(function.names[i]),
            list(psi,
                 best.params[1],
                 best.params[2]))

    models[[function.names[i]]][["PLC.predict.all"]] <-
      do.call(match.fun(function.names[i]),
              list(psi_all,
                   best.params[1],
                   best.params[2]))
  }


  return(models)
}
