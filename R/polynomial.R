#' @name polynomial
#' @title polynomial
#' @author Félicien Meunier
#' @export
#' @description Returns polynomial PLC
#' @param psi water potential
#' @param a a
#' @param b b (P50)

polynomial <- function(psi,a = 2, b = -2){
  PLC <- 100 * (1-(1 + abs(psi)^a/abs(b))^(-1))
  return(PLC)
}

#' @name polynomial.comp
#' @title polynomial.comp
#' @author Félicien Meunier
#' @export
#' @description Returns RMSE
#' @param data data.frame with water potential ("psi") and PLC ("PLC")
#' @param a a
#' @param b b (P50)
#'
polynomial.comp <- function(data,a = 2, b = -2){
  psi <- data[["psi"]]
  PLC_mes <- data[["PLC"]] ;
  PLC_mod <- polynomial(psi,a,b)
  N <- length(PLC_mes)
  RMSE <- sqrt(sum((PLC_mes-PLC_mod)^2)/(N-1))
  return(list(RMSE = RMSE,PLC_mod = PLC_mod))
}

#' @name invert.polynomial
#' @title invert.polynomial
#' @author Félicien Meunier
#' @export
#' @description Returns inverse polynomial
#' @param x PLC
#' @param a a
#' @param b b (P50)
invert.polynomial <- function(x,a = 2, b = -2){

  psi <- -abs((abs(b)*(-1 + 1/(1 - x/100))^(1/a) ))
  return(psi)
}


#' @name slope.polynomial
#' @title slope.polynomial
#' @author Félicien Meunier
#' @export
#' @description Returns slope of polynomial
#' @param x PLC value
#' @param a a
#' @param b b (P50)
slope.polynomial <- function(x,a = 2, b = -2){

  psi <- invert.polynomial(x,a,b)
  S <- (100*a*abs(psi)^(a - 1)*sign(psi))/(abs(b)*(abs(psi)^a/abs(b) + 1)^2);

  return(S)
}

