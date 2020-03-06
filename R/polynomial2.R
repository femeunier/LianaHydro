#' @name polynomial2
#' @title polynomial2
#' @author Félicien Meunier
#' @export
#' @description Returns polynomial PLC
#' @param psi water potential
#' @param a a
#' @param b b (P50)

polynomial2 <- function(psi,a = 2, b = -2){
  PLC <- 100 * (1-1/(1 + (abs(psi)/abs(b))^a))
  return(PLC)
}

#' @name polynomial2.comp
#' @title polynomial2.comp
#' @author Félicien Meunier
#' @export
#' @description Returns RMSE
#' @param data data.frame with water potential ("psi") and PLC ("PLC")
#' @param a a
#' @param b b (P50)
#'
polynomial2.comp <- function(data,a = 2, b = -2){
  psi <- data[["psi"]]
  PLC_mes <- data[["PLC"]] ;
  PLC_mod <- polynomial2(psi,a,b)
  N <- length(PLC_mes)
  RMSE <- sqrt(sum((PLC_mes-PLC_mod)^2)/(N-1))
  r2 <- summary(lm(data = data.frame(y = PLC_mod,x = PLC_mes),formula = y ~x))$r.squared
  return(list(RMSE = RMSE,PLC_mod = PLC_mod,r2 = r2))
}

#' @name invert.polynomial2
#' @title invert.polynomial2
#' @author Félicien Meunier
#' @export
#' @description Returns inverse polynomial2
#' @param x PLC
#' @param a a
#' @param b b (P50)
invert.polynomial2 <- function(x,a = 2, b = -2){

  psi <- -abs( abs(b) * ( (1-x/100)^(-1) -1 )^(1/a) )
  return(psi)
}


#' @name slope.polynomial2
#' @title slope.polynomial2
#' @author Félicien Meunier
#' @export
#' @description Returns slope of polynomial2
#' @param x PLC value
#' @param a a
#' @param b b (P50)
slope.polynomial2 <- function(x,a = 2, b = -2){

  psi <- invert.polynomial2(x,a,b)
  S <- (100*a*sign(psi)*(abs(psi)/abs(b))^(a - 1))/(abs(b)*((abs(psi)/abs(b))^a + 1)^2)

  return(S)
}

