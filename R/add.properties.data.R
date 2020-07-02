#' @name add.properties.data
#' @title add.properties.data
#' @author Félicien Meunier
#' @export
#' @description Returns data with pPs added
#' @param data data
#'
add.properties.data <- function(data,Ps = c(12,88)){

  data.mod <- data

  for (iPs in seq(1,length(Ps))){
    data.mod[[paste0("p",Ps[iPs])]] <-NA
  }

  N <- nrow(data)

  FN <- c("weibull","sigmoidal","polynomial","polynomial2","cumnorm")

  for (i in seq(1,N)){

    currentId <- data.mod$Id[i]
    currentGF <- data.mod$GrowthForm[i]
    currentP50 <- data.mod$p50[i]
    currentax <- data.mod$ax[i]
    currentfun <- data.mod$Function[i]

    if(is.finite(currentId)){
      file <- paste0("./data/",tolower(currentGF),"rawdata.csv")
      data.temp <- read.csv(file,header = TRUE) %>% mutate(psi = - abs(psi)) %>% dplyr::select(Id,psi,PLC) %>% filter(Id == currentId)

      models <- opt.data(data = data.temp, function.names = FN)
      models <- add.properties(models,x = Ps)
      best.model <- find.best.model(models)[[1]]

      for (iPs in seq(1,length(Ps))){
        data.mod[[paste0("p",Ps[iPs])]][i] <- best.model$invert[which(names(best.model$invert) == paste0("P",Ps[iPs]))]
      }

    } else if (!is.na(currentP50) & !is.na(currentax) & !is.na(currentfun) & str_length(currentfun)>0){
      params <- obs2params(ax50=-abs(currentax),P50=currentP50,funbest = currentfun)

      functionopt <- match.fun(paste0("invert.",currentfun))
      Psall <- do.call(functionopt,list(Ps,params[1],params[2]))

      for (iPs in seq(1,length(Ps))){
        data.mod[[paste0("p",Ps[iPs])]][i] <- Psall[iPs]
      }

    }
  }

  return(data.mod)
}
