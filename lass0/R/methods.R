#' Visualizing Lasso-Zero's estimates
#'
#' Plots the regression coefficients obtained by Lasso-Zero.
#' 
#' For a \code{"lass0"} object, produces boxplots of the \code{M} obtained
#' estimates for each regression coefficient and indicates the threshold level
#' \code{tau}. Coefficients whose median is larger than \code{tau} is absolute
#' value are the ones selected by Lasso-Zero. Note that if \code{lass0} was
#' called with \code{standardizeX = TRUE}, the coefficients and threshold are
#' represented on the standardized scale.
#'
#' @param x a \code{"lass0"} object
#' @param ... further arguments that can be passed to \code{plot}
#'
#' @export
#' @import graphics
#'
#' @seealso  \code{\link{lass0}} and \code{\link{qut.MC}}
#'
#' @references Descloux, P., & Sardy, S. (2018). Model selection with
#'   lasso-zero: adding straw to the haystack to better find needles. arXiv
#'   preprint arXiv:1805.05133. \url{https://arxiv.org/abs/1805.05133}

plot.lass0 <- function(x, ...) {
    # x: a fitted "lass0" object
    arg <- list(...)
    if ("xlab" %in% names(arg)) {
        xlab <- arg$xlab
    } else{
        xlab <- "feature"
    }
    if ("ylab" %in% names(arg)) {
        ylab <- arg$ylab
    } else {
        ylab <- "coefficients (on standardized scale)" 
    }
    if ("border" %in% names(arg)) {
        border <- arg$border
    } else {
        border <- (1:length(x$coefficients)) %in% x$selected + 1
    }
    
    col.ix <- (1:length(x$coefficients)) %in% x$selected + 1
    boxplot(t(x$Betas), xlab = xlab, ylab = ylab, border = border, ...)
    abline(h=c(x$tau, -x$tau), lty=2)
}








#' Print a lass0 object
#'
#' Print a summary of the Lasso-Zero estimate
#'
#' The call that produced the object \code{x} is printed, followed by the
#' estimated regression coefficients and intercept.
#'
#' @param x a \code{"lass0"} object.
#' @param ... additional print arguments.
#' 
#' @export
#'
#' @seealso  \code{\link{lass0}} and \code{\link{plot.lass0}}
#'
#' @references Descloux, P., & Sardy, S. (2018). Model selection with
#'   lasso-zero: adding straw to the haystack to better find needles. arXiv
#'   preprint arXiv:1805.05133. \url{https://arxiv.org/abs/1805.05133}
#'   
print.lass0 <- function(x, ...) {
    cat("Call:\n")
    print(x$call)
    cat("\n Coefficients:\n")
    print(x$coefficients)
    cat("\n Intercept:\n")
    print(x$intercept)
}



#' Predict method for a Lasso-Zero fit
#' 
#' Predicted values for the response given a new input matrix Xnew, based on a \code{lass0} fit.
#' 
#' @param object a \code{"lass0"} object
#' @param Xnew a new input matrix whose number of columns equals the number of coefficients returned in \code{obj}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return vector of predictions 
#' 
#' @export
#' 
#' @seealso  \code{\link{lass0}}
#'
#' @references Descloux, P., & Sardy, S. (2018). Model selection with
#'   lasso-zero: adding straw to the haystack to better find needles. arXiv
#'   preprint arXiv:1805.05133. \url{https://arxiv.org/abs/1805.05133}
#' 
predict.lass0 <- function(object, Xnew, ...){
    Xnew <- as.matrix(Xnew)
    if(ncol(Xnew) != length(object$coefficients)) stop("ncol(X) must equal length of object$coefficients")
    object$intercept + Xnew[, object$selected, drop=FALSE] %*% object$coefficients[object$selected]
}