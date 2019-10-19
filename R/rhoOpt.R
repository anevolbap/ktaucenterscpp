  #' rhoOpt function
  #'
  #' An implementation of quasi optimal rho functions following reference [1]
  #' @param x A real number.
  #' @param cc a tunning constant.
  #' @return rho(\code{x}/\code{c}\code{c})
  #' @examples
  #' val= rhoOpt(x=0.5,cc=1)
  #' @references  [1] Salibian-Barrera, M., Willems, G., & Zamar, R. (2008).
  #' The fast-tau estimator for regression. Journal of Computational and Graphical Statistics, 17(3), 659-682.
  #' @export
  #'
  #'
  #'
  rhoOpt <- function(x, cc){
    tmp <- x^2 / 2 / (3.25*cc^2)
    tmp2 <- (1.792 - 0.972 * x^2 / cc^2 + 0.432 * x^4 / cc^4 - 0.052 * x^6 / cc^6 + 0.002 * x^8 / cc^8) / 3.25
    tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
    tmp[abs(x) > 3*cc] <- 1
    tmp
  }



  #' psiOpt
  #' psi function  is the derivative of optimal rho function (see \code{\link{rhoOpt}})
  #' @param x A real number.
  #' @param cc a tunning constant.
  #' @return rho''(\code{x}/ \code{c}\code{c})
  #' @examples
  #' psiOpt(x=0.5,cc=1)
  #' @export
  psiOpt <- function(x, cc){
    tmp <- x / (3.25*cc^2)
    tmp2 <- (-1.944 * x / cc^2 + 1.728 * x^3 / cc^4 - 0.312 * x^5 / cc^6 + 0.016 * x^7 / cc^8) / 3.25
    tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
    tmp[abs(x) > 3*cc] <- 0
    tmp
  }

  #' derpsiOpt
  #' the derivative of the psi function  (see \code{\link{psiOpt}})
  #'
  #' @param x A real number.
  #' @param cc a tunning constant.
  #' @return rho'(\code{x} / \code{c}\code{c})
  #' @examples
  #' psiOpt(x=0.5,cc=1)
  #'
  #' @export
derpsiOpt<- function(x, cc){
  tmp <- rep(1,length(x)) / (3.25*cc^2)
  # derivative x^1
  # tmp2 <- (-1.944 * x / cc^2 + 1.728 * x^3 / cc^4 - 0.312 * x^5 / cc^6 + 0.016 * x^7 / cc^8) / 3.25
  # derivative of expresion above !!
  ## tmp2 <- (-1.944/ cc^2 + (3*1.728) * x^2 / cc^4 - (5*0.312) * x^4 / cc^6 + (7*0.016) * x^6 / cc^8) / 3.25
  tmp2 <- (-1.944/ cc^2 + (5.184) * x^2 / cc^4 - (1.56) * x^4 / cc^6 + ( 0.112) * x^6 / cc^8) / 3.25
  tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
  tmp[abs(x) > 3*cc] <- 0
  tmp
  }
