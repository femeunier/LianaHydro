#' @name merge.sources
#' @title merge.sources
#' @author Félicien Meunier
#' @export
#' @description Returns data with ax/p50 corrected
#' @param data data
#'

merge.sources <- function(data,tol=1){

N <- nrow(data)
Nold=2*N

while(Nold !=N){
  print("...merging...")
  merged <- c()
  data.mod <- data.frame(matrix(ncol = ncol(data)-1,nrow = 0))
  names(data.mod) <- names(data[1:5])

  for (i in seq(1,N)){
    if (!(i %in% merged )){
      dist <- sqrt(rowSums((t(matrix(rep(c(data$Long[i],data$Lat[i]),N),nrow=2))- as.matrix(data %>% select(Long,Lat)))^2))
      r <- data$r[i]
      tomerge <- which(dist < r*tol)

      data.mod <- rbind(data.mod,
                        data.frame(Long = weighted.mean(data$Long[tomerge],data$N[tomerge]),
                                   Lat = weighted.mean(data$Lat[tomerge],data$N[tomerge]),
                                   Liana = sum(data$Liana[tomerge]),
                                   Tree = sum(data$Tree[tomerge]),
                                   N = sum(data$N[tomerge])))
      merged <- c(merged, tomerge)
    }
  }

  data.mod <- data.mod %>%  mutate(r = rmin + rmax*(N-min(N))/(max(N)-min(N)))
  Nold <- N
  N <- nrow(data.mod)
  data <- map_df(data.mod, rev)

  merged <- c()
  data.mod <- data.frame(matrix(ncol = ncol(data)-1,nrow = 0))
  names(data.mod) <- names(data[1:5])

  for (i in seq(1,N)){
    if (!(i %in% merged )){
      dist <- sqrt(rowSums((t(matrix(rep(c(data$Long[i],data$Lat[i]),N),nrow=2))- as.matrix(data %>% select(Long,Lat)))^2))
      r <- data$r[i]
      tomerge <- which(dist < r*tol)

      data.mod <- rbind(data.mod,
                        data.frame(Long = weighted.mean(data$Long[tomerge],data$N[tomerge]),
                                   Lat = weighted.mean(data$Lat[tomerge],data$N[tomerge]),
                                   Liana = sum(data$Liana[tomerge]),
                                   Tree = sum(data$Tree[tomerge]),
                                   N = sum(data$N[tomerge])))
      merged <- c(merged, tomerge)
    }
  }

  data.mod <- data.mod %>%  mutate(r = rmin + rmax*(N-min(N))/(max(N)-min(N)))
  Nold <- N
  N <- nrow(data.mod)
  data <- map_df(data.mod, rev)

}

  return(data.mod)
}
