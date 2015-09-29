################################################################################
### Title: wrapperWMW.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 08/set/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 08/set/2015 - 14:13:31:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Wrapper function for Wilcoxon-Mann-Whitney test
#' 
#' Utility function performing WMW test on each OTU individually, and then correcting 
#' for multiplicity with one of the standard corrections, \code{"fdr"} is the default. 
#'
#' @title Wrapper function for Wilcoxon-Mann-Whitney test
#' 
#' @param x
#'     \code{matrix} containing the first sample
#' @param y
#'     \code{matrix} containing the second sample
#' @param ...
#'     other parameters for \link{\code{wilcox.test}}
#' @return 
#'     \code{list} with two elements of length \code{ncol(x)}, 
#'     \code{"statistic"} containing the test-statistics, and 
#'     \code{"p.value"} containing the adjusted p.values vector.
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
wrapperWMW <- function(x, y, adjMethod = "BH", ...)
{
  res <- lapply(seq_len(NCOL(x)), 
      FUN = function(i)
      {
        suppressWarnings(stats:::wilcox.test.default(
            x = x[, i], y = y[, i], ...)[c("statistic", "p.value")])
        
      })
  rawPvals <- sapply(res, elNamed, name = "p.value")
  
  list(
      "statistic" = sapply(res, elNamed, name = "statistic"),
      "p.value" = p.adjust(rawPvals, method = adjMethod)
  )
}# END: function - wrapperWMW.R

