#' @name sigmoidal
#' @title sigmoidal
#' @author Félicien Meunier
#' @export
#' @description Returns sigmoidal PLC
#' @param psi water potential
#' @param a a
#' @param b b (P50)

sigmoidal <- function(psi,a = 2, b = -2){
  PLC <-  100 * (1./(1 + exp(a*((psi) - (b)) )))
  return(PLC)
}

#' @name sigmoidal.comp
#' @title sigmoidal.comp
#' @author Félicien Meunier
#' @export
#' @description Returns RMSE
#' @param data data.frame with water potential ("psi") and PLC ("PLC")
#' @param a a
#' @param b b (P50)
#'
sigmoidal.comp <- function(data,a = 2, b = -2){
  psi <- data[["psi"]]
  PLC_mes <- data[["PLC"]] ;
  PLC_mod <- sigmoidal(psi,a,b)
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
#' @param a a
#' @param b b (P50)
invert.sigmoidal <- function(x,a = 2, b = -2){

  psi <- - abs(((log(100/x  -1 )) /a + b ))
  return(psi)
}


#' @name slope.weibull
#' @title slope.weibull
#' @author Félicien Meunier
#' @export
#' @description Returns slope of weibull
#' @param x PLC value
#' @param a a
#' @param b b (P50)
slope.sigmoidal <- function(x,a = 2, b = -2){

  psi <- invert.sigmoidal(x,a,b)
  S <- -(100*a*exp(a*(psi-b)))/(exp(a*(psi-b)) + 1)^2

  return(S)
}

