#' Variable selection for linear regression with Lasso-Zero
#'
#' Fits a (possibly high-dimensional) linear model with Lasso-Zero. Lasso-Zero
#' aggregates several estimates obtained by solving the basis pursuit problem
#' after concatenating random noise dictionaries to the input matrix. The
#' procedure is described in more details in the paper linked to in the
#' References section below.
#'
#' @param X input matrix of dimension \code{n x p}; each row is an observation
#'   vector.
#' @param y response vector of size \code{n}.
#' @param q size of noise dictionaries. A noise dictionary consists in a
#'   Gaussian matrix G of size \code{n x q} concatenated horizontally to the
#'   input matrix X. Default is \code{q = nrow(X)}.
#' @param M number of noise dictionaries used.
#' @param tau a positive threshold value. If missing, then \code{alpha} must be
#'   supplied.
#' @param alpha level of the quantile universal threshold (number between 0 and
#'   1). If missing, then \code{tau} must be supplied.
#' @param sigma standard deviation of the noise. If \code{sigma = NULL}
#'   (default) and \code{tau = NULL}, the quantile universal threshold is
#'   computed based on a pivotal statistic.
#' @param intercept whether an intercept should be fitted. If \code{TRUE}
#'   (default), \code{y} and the columns of \code{X} are mean-centered before
#'   the analysis, and the intercept is estimated by \code{mean(y) - colMeans(X)
#'   \%*\% coefficients}.
#' @param standardizeX whether the columns of \code{X} should be standardized to
#'   have unit standard deviation. Default is \code{TRUE}.
#' @param standardizeG either a positive numerical value indicating the desired
#'   Euclidean norm of all columns of the noise dictionaries, or a logical value
#'   indicating whether the columns of the noise dictionaries should be
#'   standardized to have unit standard deviation. If \code{NULL} (default),
#'   then it is set to \code{standardizeG = standardizeX}.
#' @param qut.MC.output an object of type \code{"qut.MC"} (output of
#'   \code{qut.MC} function), providing the result of Monte Carlo simulations
#'   necessary for the approximation of the Quantile Universal Threshold. By
#'   default, \code{qut.MC.output = NULL} and the \code{qut.MC} function is
#'   called unless \code{tau} is supplied.
#' @param GEVapprox whether to approximate the distribution of the null
#'   thresholding statistic by a GEV distribution (ignored if \code{tau} is
#'   supplied). Default is \code{TRUE}.
#' @param parallel if \code{TRUE}, use parallel \code{foreach} to make
#'   computations with different noise dictionaries and to perform Monte Carlo
#'   simulations for estimating the quantile universal threshold. Must register
#'   parallel beforehand, e.g. with \code{doParallel}. Default is \code{FALSE}.
#' @param soft.thresholding if \code{TRUE}, the coefficients are soft
#'   thresholded (rather than hard thresholded) at level \code{tau}. Default is
#'   \code{FALSE}.
#' @param ols whether to refit the nonzero coefficients with an ordinary least
#'   squares procedure. Default is \code{TRUE}.
#' @param ... further arguments that can be passed to \code{qut.MC}.
#'
#' @return An object of class \code{"lass0"}. It is a list containing the
#'   following components: \item{coefficients}{estimated regression
#'   coefficients.} \item{intercept}{intercept value.}
#'   \item{fitted.values}{fitted values.} \item{residuals}{residuals.}
#'   \item{selected}{set of selected features.} \item{tau}{threshold value.}
#'   \item{Betas}{matrix of size \code{p x M} containing the values of the
#'   \code{M} estimates for the regression coefficients (on the standardized
#'   scale if \code{standardizeX = TRUE}).} \item{Gammas}{matrix of size \code{q
#'   x M} containing the values of the \code{M} obtained noise coefficient
#'   vectors (on the standardized scale unless \code{standardizeG = FALSE}).}
#'   \item{madGammas}{statistics based on the noise coefficients, corresponding
#'   to the MAD of all nonzero entries in \code{Gammas}} \item{sdsX}{standard
#'   deviations of all columns of \code{X}. Can be used to transform
#'   \code{Betas} to the original scale doing \code{Betas / sdsX}.}
#'   \item{qut.MC.output}{either the list returned by \code{qut.MC}, or a character string
#'   explaining why \code{qut.MC} was not called.} \item{call}{matched call.}
#'
#' @import stats
#'
#' @examples
#' #### EXAMPLE 1: fast example with 5x10 input matrix and a small number 
#' #### (MCrep = 50) of Monte Carlo replications for computing QUT.
#' 
#' set.seed(201)
#' ## design matrix
#' n <- 5
#' p <- 10
#' X <- matrix(rnorm(n*p), n, p)
#' ## sparse vector
#' S0 <- 1:2 # support
#' beta0 <- rep(0, p)
#' beta0[S0] <- 2
#' ## response:
#' y <- X[, S0] %*% beta0[S0] + rnorm(n)
#' ## lasso-zero:
#' lass0.obj <- lass0(X, y, alpha = 0.05, MCrep = 50)
#' betahat <- lass0.obj$coefficients
#' plot(lass0.obj)
#' 
#' 
#' #### EXAMPLE 2: with 50x100 input matrix
#' 
#' \donttest{
#' set.seed(202)
#' ## design matrix
#' n <- 50
#' p <- 100
#' X <- matrix(rnorm(n*p), n, p)
#' ## sparse vector
#' S0 <- 1:3 # support
#' beta0 <- rep(0, p)
#' beta0[S0] <- 2
#' ## response:
#' y <- X[, S0] %*% beta0[S0] + rnorm(n)
#' 
#' ## 1) lasso-zero tuned by QUT with unknown noise level
#' lass0.obj1 <- lass0(X, y, alpha = 0.05)
#' betahat1 <- lass0.obj1$coefficients
#' plot(lass0.obj1)
#'
#' ## 2) lasso-zero tuned by QUT with known noise level
#' lass0.obj2 <- lass0(X, y, alpha = 0.05, sigma = 1)
#' betahat2 <- lass0.obj2$coefficients
#'
#' ## 3) lasso-zero with fixed threshold tau = 1
#' lass0.obj3 <- lass0(X, y, tau = 1)
#' betahat3 <- lass0.obj3$coefficients
#' }
#'
#' @export
#' @importFrom foreach %dopar%
#'
#' @seealso  \code{\link{qut.MC}}
#'
#' @references Descloux, P., & Sardy, S. (2018). Model selection with
#'   lasso-zero: adding straw to the haystack to better find needles. arXiv
#'   preprint arXiv:1805.05133. \url{https://arxiv.org/abs/1805.05133}

lass0 <- function(X, y, tau, alpha, q = nrow(X), M = 30, sigma = NULL,
                  intercept = TRUE, standardizeX = TRUE,  standardizeG = NULL, 
                  qut.MC.output = NULL, GEVapprox = TRUE,
                  parallel = FALSE, soft.thresholding = FALSE, ols = TRUE, ...) {
    
    ## checking for some errors in the arguments & warnings:
    if (missing(tau)) {
        if (missing(alpha)) stop("tau and alpha are both missing; one of them must be provided")
        if (!is.numeric(alpha) | alpha > 1 | alpha < 0) stop("alpha must be a number between 0 and 1")
        if (!is.null(qut.MC.output)) {
            if (class(qut.MC.output) != "qut.MC") stop("qut.MC.output must be an object of class 'qut.MC'")
            lass0settings <- list(q = q, M = M, sigma = sigma, intercept = intercept,
                                  standardizeX = standardizeX, standardizeG = standardizeG)
            if (!identical(lass0settings, qut.MC.output$lass0settings)) {
                warning ("Some of the arguments q, M, sigma, intercept, standardizeX and standardizeG 
                         are not identical to the ones passed to qut.MC for computation of qut.MC.output
                         --> this might lead to an unappropriate threshold tau and affect your results!")
            }
        }
    } else {
        if (tau < 0) stop("tau must be positive")
        if (!missing(alpha)) stop("two many arguments: tau and alpha cannot be both supplied")
        if (!is.null(qut.MC.output)) warning("value for tau is supplied --> qut.MC.output is ignored!")
        if (!is.null(sigma)) warning("value for tau is supplied --> sigma is ignored!")
        qut.MC.output <- "value for tau is supplied--> qut.MC was not called" 
    }
    if (M < 1) {
        stop("M must be a number >= 1")
    } else if (M > 1 & q == 0) {
        warning("q = 0 --> no noise dictionary is used, therefore M is set to 1")
        M <- 1
    } 
    
    #
    X <- as.matrix(X)
    p <- ncol(X)
    n <- nrow(X)
    if (is.null(standardizeG)) standardizeG <- standardizeX
    if (intercept) {
        meany <- mean(y)
        y <- y - meany
        meansX <- colMeans(X)
        X <- t(t(X) - colMeans(X))
    } else {
        meany <- 0
    }
    if (standardizeX) {
        sdsX <- apply(X, 2, sd)
        X <- t(t(X) / sdsX)
    } else {
        sdsX <- rep(1, p)
    }
    
    ## extended Basis Pursuit with noise dictionaries:
    extBP <- MextBP(X, y, q = q, M = M, meancenterG = intercept, standardizeG, 
                     parallel = parallel, returnGammas = TRUE)
    Betas <- extBP$Betas
    Gammas <- extBP$Gammas
    madGammas <- extBP$madGammas
    betahat <- extBP$medBeta
    features <- colnames(X)
    if (is.null(features)) features <- paste0("X", 1:ncol(X))
    if(is.vector(Betas)) {
        names(Betas) <- features
    } else {
        rownames(Betas) <- features
    }
    names(betahat) <- features
    
    ## thresholding:
    if (missing(tau)) {
        if (is.null(qut.MC.output)) {
            print("MC simulations for QUT estimation")
            qut.MC.output <- qut.MC(X = X, q = q, M = M, alpha = alpha, sigma = sigma,
                          intercept = intercept, standardizeX = FALSE, 
                          standardizeG = standardizeG, GEVapprox = GEVapprox,
                          parallel = parallel, ...)
            upperQuant <- qut.MC.output$upperQuant
            if (is.null(sigma)) {
                tau <- madGammas * upperQuant
            } else {
                tau <- upperQuant
            }
        } else {
            if (GEVapprox) {
                GEVpar <- qut.MC.output$GEVpar
                if(is.character(GEVpar)) {
                    warning("qut.MC.output contains no GEV parameters --> they are estimated from qut.MC.output$allMC")
                    print("fitting GEV distribution")
                    GEVfit <- ismev::gev.fit(qut.MC.output$allMC)
                    GEVpar <- GEVfit$mle
                }
                if (GEVpar[3] == 0) {
                    upperQuant <- GEVpar[1] - GEVpar[2] * log(-log(1-alpha))
                } else {
                    upperQuant <- GEVpar[1] + GEVpar[2] / GEVpar[3] * (1/(-log(1-alpha))^GEVpar[3] - 1)
                }
            } else {
                upperQuant <- quantile(qut.MC.output$allMC, 1-alpha)
            }
            if (is.null(sigma)) {
                tau <- madGammas * upperQuant
            } else {
                tau <- upperQuant
            }
        }
    }
    
    betahat[abs(betahat) <= tau] <- 0
    if (soft.thresholding) {
        ix <- which(abs(betahat) > tau)
        betahat[ix] = betahat[ix] - sign(betahat[ix]) * tau
    }
    selected <- which(betahat != 0)
    if (length(selected) > 0 & ols) {
        betahat[selected] <- lm(y ~ X[, selected, drop = FALSE] - 1)$coefficients
    }
    
    ## transforming to original scale, get coefficients, fitted values etc.
    coefficients <- betahat / sdsX
    if (intercept) {
        intercept <- as.numeric(meany - t(meansX) %*% coefficients)
    } else {
        intercept <- 0
    }
    fitted.values <- meany + X[, selected, drop = FALSE] %*% betahat[selected]
    residuals <- (y + meany) - fitted.values
    
    ## "lass0" object to return:
    out <- list()
    out$coefficients <- coefficients
    out$intercept <- intercept
    out$fitted.values <- fitted.values
    out$residuals <- residuals
    out$selected <- selected
    out$tau <- tau
    out$Betas <- Betas
    out$Gammas <- Gammas
    out$madGammas <- madGammas
    out$sdsX <- sdsX
    out$qut.MC.output <- qut.MC.output
    out$call <- match.call()
    class(out) <- "lass0"
    out
}