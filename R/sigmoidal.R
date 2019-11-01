#' @name sigmoidal
#' @title sigmoidal
#' @author Félicien Meunier
#' @export
#' @description Returns sigmoidal PLC
#' @param psi water potential
#' @param param vector of a and b (P50)

sigmoidal <- function(psi,param){
  a <- param[1]
  b <- param[2]
  PLC <-  100 * (1./(1 + exp(a*((psi) - (b)) )))
  return(PLC)
}

#' @name sigmoidal.comp
#' @title sigmoidal.comp
#' @author Félicien Meunier
#' @export
#' @description Returns RMSE
#' @param data data.frame with water potential ("psi") and PLC ("PLC")
#' @param param vector of a and b
#'
sigmoidal.comp <- function(data,param){
  psi <- data[["psi"]]
  PLC_mes <- data[["PLC"]] ;
  PLC_mod <- sigmoidal(psi,param)
  N <- length(PLC_mes)
  RMSE <- sqrt(sum((PLC_mes-PLC_mod)^2)/(N-1))
  return(list(RMSE = RMSE,PLC_mod = PLC_mod))
}

#' @name invert.sigmoidal
#' @title invert.sigmoidal
#' @author Félicien Meunier
#' @export
#' @description Returns inverse sigmoidal
#' @param x PLC
#' @param param vector of a and b
#'
invert.sigmoidal <- function(x,param){
  a <- param[1]
  b <- param[2]

  psi <- - abs(((log(100/x  -1 )) /a + b ))
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
slope.sigmoidal <- function(x,param){

  a <- param[1]
  b <- param[2]

  psi <- invert.sigmoidal(x,param)
  S <- -(100*a*exp(-a*(b - psi)))/(exp(-a*(b - psi)) + 1)^2

  return(S)
}

