################################################################################
### Title: internals.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 04/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 04/mar/2015 - 16:21:48:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################

###### my version of Weir computation (Method of Moments)
weirMoM4Wald <- function (data, se = FALSE)
{
  if (missing(data)) 
    stop("data missing.")
  K <- ncol(data)
  J <- nrow(data)
  
  totalsample <- sum(data, na.rm = TRUE)
  MoM <- colSums(data, na.rm = TRUE)/totalsample
  Sn <- rowSums(data, na.rm = TRUE)
  auxData <- data / Sn
  
  MSP <- sum(rowSums(
          (auxData - matrix(rep(MoM, J), nrow = J, ncol = K, byrow = TRUE))^2,
      na.rm = TRUE) * Sn) / (J - 1)
  MSG <- sum(rowSums(auxData * (1 - auxData), na.rm = TRUE) * Sn) / (totalsample - J)
  
  nc <- (sum(Sn) - sum(Sn^2)/sum(Sn)) / (J - 1)
  MoM.wh <- (MSP - MSG)/(MSP + (nc - 1) * MSG)
  
  if (se)
  {
    std.er <- sqrt(2 * (1 - MoM.wh)^2 / (J - 1) * 
            ((1 + (nc - 1) * MoM.wh)/nc)^2)
    return(list(theta = MoM.wh, se = std.er))
  } else
  {
    return(MoM.wh)
  }
}# END - myWeirMom


### pi vector estimation with Method of Moments
piMoM4Wald <- function(data)
{
  totalReads <- sum(data, na.rm = TRUE)
  piMom <- colSums(data, na.rm = TRUE)/totalReads
  zeroInds <- abs(piMom) < .Machine$double.eps
  r <- sum(zeroInds)
  rr <- length(piMom) - r
  piMom[!zeroInds] <- piMom[which(piMom != 0)] - r/(rr * 2 * (totalReads + 1))
  piMom[zeroInds] <- 1/(2 * (totalReads + 1))
  
  return(piMom)
}


