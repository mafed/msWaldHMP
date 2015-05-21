################################################################################
### Title: msWald.R                                              ###
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
###     > 0.1 - 04/mar/2015 - 16:27:52:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Wrapper for msWaldStat
#' 
#' @title Wrapper for msWaldStat
#' 
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
#' @param wmnTest 
#'   \code{logical} default to \code{FALSE}, if \code{TRUE} performs the two-samples 
#'   Wilcoxon-Mann-Whitney test one OTU at-a-time and then corrects for multiplicity.
#'   {\bfseries Note:} it works only with case (i), thus without subsets and 
#'   without strata.
#' @param adjMethod
#'   \code{character} method to adjust for multiplicity in case \code{wmnTest} is 
#'   \code{TRUE}
#' @return 
#'     what the function returns
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
### nReads is a list of number of reads, i-th position of the list contains 
### n_i elements: number of reads of that observation/sample (library size)
msWald <- function(nReads, alphaDM, thetaDM, wmnTest = FALSE, adjMethod = "fdr")
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
  
  ### checks of thetaVec
  thetaDM <- as.matrix(thetaDM)
  if (ncol(thetaDM) < nStrata)
  {
    thetaDM <- matrix(thetaDM, nrow = length(thetaDM), ncol = nStrata)
#    message("  theta parameters has been recycled")
  } else {}
  rownames(thetaDM) <- paste0("group", seq_len(nGroups))
  colnames(thetaDM) <- names(alphaDM)
  
  
  sampleSizes <- sapply(nReads, FUN = function(x) sapply(x, FUN = length))
  rownames(sampleSizes) <- paste0("group", seq_len(nGroups))
  colnames(sampleSizes) <- names(nReads)
  maxSampleSize <- apply(sampleSizes, MARGIN = 1, max)
  
  dmDataList <- vector("list", nGroups)
  names(dmDataList) <- paste0("group", seq_len(nGroups))
  
  
  piMomAux <- array(NA, dim = c(dim(el(alphaDM)), nSubsets), 
          dimnames = list(
              rownames(el(alphaDM)), colnames(el(alphaDM)), 
              names(nReads))
      )
  
  thetaMomAux <- matrix(NA, nrow = nGroups, ncol = nSubsets, 
      dimnames = list(
          paste0("group", seq_len(nGroups)), 
          names(nReads)
      ))
  
  res <- matrix(NA, nrow = nStrata, ncol = nSubsets, 
      dimnames = list(names(alphaDM), names(nReads)))
  
  
  aux <- nReads[apply(sampleSizes, 1, which.max)]
  nReadsMax <- lapply(seq_along(aux), 
      FUN = function(i, dat) dat[[c(i, i)]], dat = aux)
  
  
  ### main loop over strata (enterotypes), usually equal to 1 or 3
  for (stRun in seq_len(nStrata))
  {
    ## loop over groups (can be more than 2)
    for (grRun in seq_len(nGroups))
    {
      dmDataList[[grRun]] <- t(sapply(seq_len(maxSampleSize[grRun]), 
              FUN = function(i, nRds, alMat)
              {
                rmultinom(n = 1, size = nRds[i], prob = rDirichlet(n = 1, alpha = alMat))
              }, 
              nRds = nReadsMax[[grRun]], 
              alMat = alphaDM[[stRun]][, grRun] * 
                  (1 - thetaDM[grRun, stRun])/thetaDM[grRun, stRun]
          ))
      
      ## populate *theta* and *pi* for the Wald statistic
      for (subRun in seq_len(nSubsets))
      {
        subsetData <- dmDataList[[grRun]][seq_len(sampleSizes[grRun, subRun]), ]
        thetaMomAux[grRun, subRun] <- msWaldHMP:::weirMoM4Wald(subsetData)
        piMomAux[, grRun, subRun] <- msWaldHMP:::piMoM4Wald(subsetData)
      }# END - for: subRun
    }# END - for: grRun
    
    ## compute Wald statistics for each subset in the stratum
    for (subRun in seq_len(nSubsets))
    {
      res[stRun, subRun] <- msWaldHMP:::msWaldStat(
          nReads[[subRun]], piMat = piMomAux[, , subRun], thetaVec = thetaMomAux[, subRun])
    }# END - for: subRun
  }# END - for: stRun 
  
  
  
  if (length(el(nReads)) == 2L && wmnTest)
  {
    nReads <- nReads[[1L]]
    alphaDM <- alphaDM[[1L]]
    
    pVals <- sapply(seq_len(nrow(alphaDM)), FUN = function(i, dat1, dat2)
        {
          wilcox.test(x = dat1[, i], y = dat2[, i], exact = FALSE)$p.value 
        }, 
        dat1 = dmDataList[[1L]], dat2 = dmDataList[[2L]]
    )
    
    pValsAdj <- p.adjust(pVals, method = adjMethod)
    res <- list(out, cbind("pVals" = pVals, "pValsAdj" = pValsAdj))
  } else {}
  
  res
  
}# END - msWald

