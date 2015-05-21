################################################################################
### Title: msWaldStat.R                                              ###
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
###     > 0.1 - 04/mar/2015 - 16:22:38:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Multi-sample Wald-type test statistics, basically Xmcupo.sevsample()
#' 
#' @title Multi-Sample Wald Statistic
#' 
#' @param nReads 
#'   list of number of reads (library size) for each group, 
#'   i-th element of the list contains \eqn{n_i} elements
#' @param piMat 
#'   matrix K x G, K = number of OTUs, G = number of groups, contains probability 
#'   vectors/R.A.D.s
#' @param thetaVec 
#'   overdispersion parameter vector, one for each group
#' @return 
#'     \code{numeric} Wald-type test statistic
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
msWaldStat <- function(nReads, piMat, thetaVec)
{
  nGroups <- length(nReads)
  
  N_1Cj <- rep.int(NA, nGroups)
  for (grRun in seq_len(nGroups))
  {
    N_1Cj[grRun] <- (thetaVec[grRun] * (sum(nReads[[grRun]]^2) - sum(nReads[[grRun]])) + 
          sum(nReads[[grRun]])) / sum(nReads[[grRun]])^2
  }
   
  N_1CjInv <- 1 / N_1Cj
  
  pi0Est <- drop((piMat %*% N_1CjInv) / sum(N_1CjInv))
  
  Xmcupo <- sum(N_1CjInv * colSums(
      (piMat - matrix(pi0Est, nrow = length(pi0Est), ncol = nGroups))^2 / drop(pi0Est)
  ))
  
  Xmcupo
}# END - msWaldStat
