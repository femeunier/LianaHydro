#' @name weibull
#' @title weibull
#' @author Félicien Meunier
#' @export
#' @description Returns Weibull PLC
#' @param psi water potential
#' @param lambda lambda
#' @param k k
weibull <- function(psi,lambda = 2, k = 2){

  PLC <- 100*(1-exp(-((abs(psi)/lambda)^k)))
  return(PLC)
}

#' @name weibul.comp
#' @title weibul.comp
#' @author Félicien Meunier
#' @export
#' @description Returns RMSE
#' @param data data.frame with water potential ("psi") and PLC ("PLC")
#' @param lambda lambda
#' @param k k
weibull.comp <- function(data,lambda = 2, k = 2){
  psi <- data[["psi"]]
  PLC_mes <- data[["PLC"]] ;
  PLC_mod <- weibull(psi,lambda,k)
  N <- length(PLC_mes)
  RMSE <- sqrt(sum((PLC_mes-PLC_mod)^2)/(N-1))
  r2 <- summary(lm(data = data.frame(y = PLC_mod,x = PLC_mes),formula = y ~x))$r.squared
  return(list(RMSE = RMSE,PLC_mod = PLC_mod,r2 = r2))
}

#' @name invert.weibull
#' @title invert.weibull
#' @author Félicien Meunier
#' @export
#' @description Returns inverse weibull
#' @param x PLC
#' @param lambda lambda
#' @param k k
invert.weibull <- function(x,lambda = 2, k = 2){

  psi <- -abs((-log(1 - x/100))^(1/k)*lambda)
  return(psi)
}


#' @name slope.weibull
#' @title slope.weibull
#' @author Félicien Meunier
#' @export
#' @description Returns slope of weibull
#' @param x PLC value
#' @param lambda lambda
#' @param k k
slope.weibull <- function(x,lambda = 2, k = 2){


  psi <- invert.weibull(x,lambda,k)
  S <- (100*k*exp(-(abs(psi)/lambda)^k)*sign(psi)*(abs(psi)/lambda)^(k - 1))/lambda

  return(S)
}

