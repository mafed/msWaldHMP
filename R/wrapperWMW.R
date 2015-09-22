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
#'     \code{numeric} vector of length \code{ncol(x)} containing adjusted p.values
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
wrapperWMW <- function(x, y, adjMethod = "BH", ...)
{
  rawPvals <- sapply(seq_len(NCOL(x)), 
      FUN = function(i)
      {
        suppressWarnings(stats:::wilcox.test.default(
            x = x[, i], y = y[, i], ...)$p.value)
        
      })
  p.adjust(rawPvals, method = adjMethod)
}# END: function - wrapperWMW.R

