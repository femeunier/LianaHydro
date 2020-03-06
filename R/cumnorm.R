#' @name cumnorm
#' @title cumnorm
#' @author Félicien Meunier
#' @export
#' @description Returns cumnorm PLC
#' @param psi water potential
#' @param m mean (P50)
#' @param std standard deviation

cumnorm <- function(psi,m = -2, std = 1){
  # std <- abs(std)
  PLC <- 100*(1-pnorm((psi),m,std))
  return(PLC)
}

#' @name cumnorm.comp
#' @title cumnorm.comp
#' @author Félicien Meunier
#' @export
#' @description Returns RMSE
#' @param data data.frame with water potential ("psi") and PLC ("PLC")
#' @param m mean (P50)
#' @param std standard deviation
#'
cumnorm.comp <- function(data,m = -2, std = 1){
  psi <- data[["psi"]]
  PLC_mes <- data[["PLC"]] ;
  PLC_mod <- cumnorm(psi,m,std)
  N <- length(PLC_mes)
  RMSE <- sqrt(sum((PLC_mes-PLC_mod)^2)/(N-1))
  r2 <- summary(lm(data = data.frame(y = PLC_mod,x = PLC_mes),formula = y ~x))$r.squared
  return(list(RMSE = RMSE,PLC_mod = PLC_mod,r2 = r2))
}

#' @name invert.cumnorm
#' @title invert.cumnorm
#' @author Félicien Meunier
#' @export
#' @description Returns inverse polynomial
#' @param x PLC
#' @param m mean (P50)
#' @param std standard deviation
invert.cumnorm <- function(x,m = -2, std = 1){

  # std <- abs(std)
  psi <- qnorm((1-x/100),m,std)
  return(psi)
}


#' @name slope.cumnorm
#' @title slope.cumnorm
#' @author Félicien Meunier
#' @export
#' @description Returns slope of cumnorm
#' @param x PLC value
#' @param m mean (P50)
#' @param std standard deviation
slope.cumnorm <- function(x,m = -2, std = 1){

  # std <- abs(std)
  psi <- invert.cumnorm(x,m,std)
  S <- 100*dnorm(psi, m, std)

  return(S)
}

