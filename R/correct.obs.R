#' @name correct.obs
#' @title correct.obs
#' @author Félicien Meunier
#' @export
#' @description Returns data with ax/p50 corrected
#' @param data data
#'
correct.obs <- function(data){

  data.mod <- data

  N <- nrow(data)

  FN <- c("weibull","sigmoidal","polynomial","polynomial2","cumnorm")

  for (i in seq(1,N)){

    currentId <- data.mod$Id[i]
    currentGF <- data.mod$GrowthForm[i]

    if(is.finite(currentId)){
      file <- paste0("./data/",tolower(currentGF),"rawdata.csv")
      data.temp <- read.csv(file,header = TRUE) %>% mutate(psi = - abs(psi)) %>% dplyr::select(Id,psi,PLC) %>% filter(Id == currentId)

      models <- opt.data(data = data.temp, function.names = FN)
      models <- add.properties(models,x = 50)
      best.model <- find.best.model(models)[[1]]

      data.mod$p50[i] <- best.model$invert[which(names(best.model$invert) == "P50")]
      data.mod$ax[i] <- abs(best.model$slopes[which(names(best.model$slopes) == "ax50")])
    }
  }

  return(data.mod)
}
