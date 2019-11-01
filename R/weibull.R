#' @name weibull
#' @title weibull
#' @author Félicien Meunier
#' @export
#' @description Returns Weibull PLC
#' @param psi water potential
#' @param param vector of lambda and k

weibull <- function(psi,param){
  lambda <- param[1]
  k <- param[2]
  PLC <- 100*(1-exp(-((abs(psi)/lambda)^k)))
  return(PLC)
}

#' @name weibul.comp
#' @title weibul.comp
#' @author Félicien Meunier
#' @export
#' @description Returns RMSE
#' @param data data.frame with water potential ("psi") and PLC ("PLC")
#' @param param vector of lambda and k
#'
weibul.comp <- function(data,param){
  psi <- data[["psi"]]
  PLC_mes <- data[["PLC"]] ;
  PLC_mod <- weibull(psi,param)
  N <- length(PLC_mes)
  RMSE <- sqrt(sum((PLC_mes-PLC_mod)^2)/(N-1))
  return(list(RMSE = RMSE,PLC_mod = PLC_mod))
}

#' @name invert.weibull
#' @title invert.weibull
#' @author Félicien Meunier
#' @export
#' @description Returns inverse weibull
#' @param x PLC
#' @param param vector of lambda and k
#'
invert.weibull <- function(x,param){
  lambda <- param[1]
  k <- param[2]

  psi <- -abs((-log(1 - x/100))^(1/k)*lambda)
  return(psi)
}


#' @name slope.weibull
#' @title slope.weibull
#' @author Félicien Meunier
#' @export
#' @description Returns slope of weibull
#' @param x PLC value
#' @param param vector of lambda and k
#'
slope.weibull <- function(x,param){

  lambda <- param[1]
  k <- param[2]

  psi <- invert.weibull(x,param)
  S <- (100*k*exp(-(abs(psi)/lambda)^k)*sign(psi)*(abs(psi)/lambda)^(k - 1))/lambda

  return(S)
}

