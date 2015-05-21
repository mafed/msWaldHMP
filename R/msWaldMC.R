################################################################################
### Title: msWaldMC.R                                              ###
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
###     > 0.1 - 04/mar/2015 - 16:31:08:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Wrapper function for MonteCarlo simulations
#' 
#' @title MonteCarlo Simulations with Multi-Sample Wald Statistic
#' 
#' @param MC 
#'   number of MonteCarlo iterations
#' @param nReads 
#'   either (i) a \code{list} or (ii) a \code{list} of \code{list}s; more in detail:
#'   \itemize{
#'     \item{}{
#'     in case (i) \code{nReads} has to be a list of number of reads (library size) for
#'     each group, {\em i.e.} the i-th element of the list contains \eqn{n_i} elements;
#'     }
#'     \item{}{
#'     in case (ii) \code{nReads} must be a \code{list} where each element contains a 
#'     \code{list} with the same structure of the previous point, {\em i.e.} each element
#'     element is used to perform a separate analysis on the same data by subsetting the
#'     samples. It useful in case of simulations when several sample sizes are to be 
#'     tried.
#'     }
#'   }
#'   See also the {\em details} section.
#' @param alphaDM 
#'   \eqn{\alpha}{alpha} parameter for the Dirichlet-Multinomial distribution.
#'   Either a \code{matrix} of dimensions \eqn{K \times G}{K x G}, K = number of OTUs, 
#'   G = number of groups, containing the parameter vectors for the Dirichlet 
#'   distribution or a \code{list} where each element contains a \eqn{K \times G}{K x G}
#'   \code{matrix} in case the simulations need to be stratified, {\em e.g.} by 
#'   enterotype. In the latter case simulations are performed separately 
#'   stratum-by-stratum; individual and global rejection rates (power) are given as 
#'   output (see also details).
#' does NOT need to be stratified (subsets of samples, 
#'     {\em e.g.} in case of enterotypes)
#' @param thetaDM 
#'   \eqn{\theta}{theta} overdispersion parameter for the Dirichlet-Multinomial 
#'   distribution. Either (i) a \code{numeric} vector of length equal 
#'   to the number of groups under test or (ii) a \code{matrix} where each column 
#'   corresponds to a different stratum. It is recycled when possible (same values for 
#'   each stratum).
#' @param sigLev 
#'   significance level (alpha, or type-I error)
#' @param avgRej 
#'   if FALSE it returns the number of rejections instead of the proportion 
#'   (among MC iterations)
#' @param ...
#'   see \code{\link{msWald}} for details
#' @return 
#'     either number of rejections or the mean
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
msWaldMC <- function(MC = 100, nReads, alphaDM, thetaDM, sigLev = .05, avgRej = TRUE,
    wmnTest = FALSE, adjMethod = "fdr")
{
  ### check class and names of *nReads*
  if (!is.list(el(nReads)))
  {
    nReads <- list("totSet" = nReads)
  } else
  {
    if (is.null(names(nReads)))
      names(nReads) <- paste0("subset", seq_along(nReads))
  }
  
  ### check class and names of *alphaList*
  if (!is.list(alphaDM))
  {
    alphaDM <- list(alphaDM)
  } else
  {
    if (is.null(names(alphaDM)))
      names(alphaDM) <- paste0("stratum", seq_along(alphaDM))
  }
  
  nGroups <- length(el(nReads))
  nStrata <- length(alphaDM)
  nSubsets <- length(nReads)
  nOtus <- sapply(alphaDM, FUN = nrow)
  
  ### quantiles of the reference distribution
  ## Degrees of Freedom
  DoFs <- (nGroups - 1) * (nOtus - 1)
  qAlpha <- qchisq(p = 1 - sigLev/nStrata, df = DoFs, ncp = 0, lower.tail = TRUE)
  qAlphaGlob <- qchisq(p = 1 - sigLev, df = sum(DoFs), ncp = 0, lower.tail = TRUE)
  
  ### MonteCarlo simulations
  tmp <- lapply(seq_len(MC), 
      FUN = function(x, wmn, adj)
      {
        msWald(nReads, alphaDM, thetaDM, wmnTest = wmn, adjMethod = adj)
      }, wmn = wmnTest, adj = adjMethod)
  
  res <- array(unlist(tmp), dim = c(nStrata, nSubsets, MC), 
      dimnames = list(
          rownames(el(tmp)), colnames(el(tmp)), paste0("MC", seq_len(MC))
      ))
  
  
  ### TODO: use *pchisq* plus multiplicity correction instead of just checking rejection
#  resAdj <- pchisq(res, df = DoFs[1L], lower.tail = FALSE)
#  resAdj <- apply(resAdj, MARGIN = 2L:3, fun = p.adjust, method = adjMethod)
#  rejRes <- rowSums(resAdj <= sigLev, na.rm = TRUE, dims = 2)
  
  
  ### number of rejections, note that the following command is equivalent to:
  ### rejRes <- apply(res > qAlpha, MARGIN = 1L:2, FUN = sum) 
  rejRes <- rowSums(res > qAlpha, na.rm = TRUE, dims = 2)
  
  if (avgRej)
  {
    rejRes <- rejRes / MC
  } else {}# END - ifelse: avgRej
  
  
  ### if multilple strata, sum up together
  if (nStrata > 1L)
  {
##   TODO: sums wald stat and compare with different quatile of chisquare
#    rejResGlob <- rowSums(apply(res > qAlpha, MARGIN = c(2L, 3), 
#        FUN = function(x) !all(!x)))
    rejResGlob <- rowSums(colSums(res, na.rm = TRUE) > qAlphaGlob)
  
    if (avgRej)
    {
      rejResGlob <- rejResGlob / MC
    }
    
    ## bind together the two
    rejRes <- rbind(rejRes, "Global" = rejResGlob)
  } else {}
  
  
  
  
  
  
  if (!is.null(dim(rejRes)) && wmnTest)
  {
    waldRes <- unlist(res[2 * seq_len(MC) - 1L]) > qAlpha
    
    wmnRes <- sapply(res[2 * seq_len(MC)], function(x, al)
        {
          colSums(x <= al, na.rm = TRUE) > 0
        }, al = sigLev)
    
    if(avgRej)
    {
      c("wald" = mean(waldRes, na.rm = TRUE), "WMN" = rowMeans(wmnRes, na.rm = TRUE))
    } else
    {
      c("wald" = sum(waldRes, na.rm = TRUE), "WMN" = rowSums(wmnRes, na.rm = TRUE))
    }# END - ifelse: avgRej
  }# END - ifelse: wmnTest
  
  
  rejRes
  
}

#set.seed(12345)
#ha <- msWaldMC(MC = 500, nReads = nReads, alphaDM = alphaDM, thetaDM = thetaDM)
#set.seed(12345)
#ha2 <- msWaldMC(MC = 500, nReads = nReads[[3L]], alphaDM = alphaDM, thetaDM = thetaDM)
#set.seed(12345)
#ha3 <- msWaldMC(MC = 500, nReads = nReads, alphaDM = alphaDM$ent1, thetaDM = thetaDM[, 1])
#
#ha
#ha2
#ha3


