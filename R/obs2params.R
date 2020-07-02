#' @name obs2params
#' @title obs2params
#' @author Félicien Meunier
#' @export
#' @description Returns function parameters
#' @param ax50 ax50
#' @param P50 P50
#' @param funbest funbest
#'
obs2params <- function(ax50,P50,funbest){

  default <- default.model.params(funbest)

  obs <- c(P50,ax50)

  retrieve_fun.param <- function(params,obs,funbest){

    # print(params)
    functionopt <- match.fun(funbest)
    functionSopt <- match.fun(paste0("slope.",funbest))
    functionIopt <- match.fun(paste0("invert.",funbest))

    P50.mod <-  do.call(functionIopt,list(50,params[1],params[2]))
    ax50.mod <-  do.call(functionSopt,list(50,params[1],params[2]))

    P50  <- obs[1]
    ax50 <- obs[2]

    diff <- sum(abs(c((ax50-ax50.mod)/ax50,(P50-P50.mod)/P50)))

    return(abs(diff))
  }


  if(sign(default[[2]]) < 0){
    xmin <- optim(c(default[[1]],default[[2]]),
                  retrieve_fun.param,
                  obs=obs, funbest=funbest,
                  lower = c(0.001,Inf*sign(default[[2]])),
                  method = c("L-BFGS-B"))
  } else{
    xmin <- optim(c(default[[1]],default[[2]]),
                  retrieve_fun.param,
                  obs=obs, funbest=funbest,
                  lower = c(0.001,0.001),
                  upper = c(Inf,Inf*sign(default[[2]])),
                  method = c("L-BFGS-B"))
  }

  params <- xmin$par
  return(params)
}
