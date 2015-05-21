################################################################################
### Title: rDirichlet.R                                              ###
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
###     > 0.1 - 04/mar/2015 - 16:20:08:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Draws from Dirichlet distribution
#' 
#' @title Random Dirichlet Draws
#' 
#' @param n
#'     number of draws
#' @param alpha
#'     parameter vector
#' @return 
#'     what the function returns
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 

rDirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

