#' Monte Carlo simulation for estimating the Quantile Universal Threshold
#'
#' Performs a Monte Carlo simulation to estimate the distribution of the null
#' thresholding statistic required for computation of the quantile universal
#' threshold, and computes its upper alpha-quantile if alpha is provided.
#'
#' If the noise level \code{sigma} is known, the statistic of interest is simply
#' the sup-norm of the Lasso-Zero coefficients obtained under the null
#' hypothesis (i.e. when all coefficients all zero) when the threshold
#' \code{tau} is set to 0, and its upper alpha-quantile is the quantile
#' universal threshold. If \code{sigma = NULL} (sigma unknown) a pivotized
#' statistic is used, which is obtained by dividing the statistic described
#' above by the MAD of all nonzero noise coefficients obtained by Lasso-Zero.
#'
#' @param X input matrix of dimension \code{n x p}, previously mean-centered and
#'   standardized if necessary!
#' @param alpha level of the quantile universal threshold. By default
#'   \code{alpha = NULL} and no quantile is returned.
#' @param sigma standard deviation of the noise. If \code{sigma = NULL}
#'   (default), the statistic of interest if pivotized.
#' @param q size of noise dictionaries. A noise dictionary consists in a
#'   Gaussian matrix G of size \code{n x q} concatenated horizontally to the
#'   input matrix X. Default is \code{q = nrow(X)}.
#' @param M number of noise dictionaries used.
#' @param intercept if \code{TRUE} (default), the columns of \code{X} are
#'   mean-centered before the analysis.
#' @param standardizeX whether the columns of \code{X} should be standardized to
#'   have unit standard deviation. Default is \code{TRUE}.
#' @param standardizeG either a positive numerical value indicating the desired
#'   Euclidean norm of all columns of the noise dictionaries, or a logical value
#'   indicating whether the columns of the noise dictionaries should be
#'   standardized to have unit standard deviation. If \code{NULL} (default),
#'   then it is set to \code{standardizeG = standardizeX}.
#' @param MCrep number of Monte Carlo replications. Default is \code{MCrep =
#'   100.}
#' @param GEVapprox whether to approximate the distribution of the null
#'   thresholding statistic by a GEV distribution. If \code{TRUE}, the maximum
#'   likelihood estimates of the GEV parameters are computed on the Monte Carlo
#'   sample. Default if \code{TRUE}.
#' @param parallel if \code{TRUE}, use parallel \code{foreach} to perform the
#'   Monte Carlo simulation. Must register parallel beforehand, e.g. with
#'   \code{doParallel}. Default is \code{FALSE}.
#'
#' @return An object of class \code{"qut.MC"}, which is a list with the
#'   following components: \item{allMC}{all \code{MCrep} realizations of the
#'   null thresholding statistic of interest (pivotized if \code{sigma =
#'   NULL}).} \item{GEVpar}{MLE estimates of the GEV distribution parameters
#'   (\code{NULL} if \code{GEVapprox} was set to \code{FALSE}).}
#'   \item{upperQuant}{upper alpha-quantile of the null thresholding statistis
#'   (either the empirical quantile, or the quantile of the fitted GEV
#'   distribution).} \item{call}{matched call.} \item{lass0settings}{a list
#'   containing the chosen settings for the computation of lasso-zero: \code{q},
#'   \code{M}, \code{sigma}, \code{intercept}, \code{standardizeX} and
#'   \code{standardizeG}. When an object of type \code{"qut.MC"} is supplied to
#'   the \code{lass0} function, a warning message appears if the corresponding
#'   arguments passed to \code{lass0} are different.}
#'
#' @export
#' @import stats
#' @importFrom foreach %dopar%
#'
#' @examples
#' ### Fast toy example with 5x10 input matrix and a small number (MCrep = 50)
#' ### of Monte Carlo replications.
#' ### Illustrates how to tune Lasso-Zero with QUT for the same input matrix but
#' ### different responses and/or different alpha values, without calling
#' ### qut.MC several times:
#'
#' ### (for faster computation when X and MCrep are larger: register a parallel 
#' ### backend and choose parallel = TRUE when calling lass0 and qut.MC functions.)
#'
#' set.seed(3)
#'
#' ## input matrix:
#' n <- 5
#' p <- 10
#' X <- matrix(rnorm(n*p), n, p)
#'
#' ## two sparse vectors and corresponding responses:
#' S1 <- 1:2 # first support
#' beta1 <- numeric(p)
#' beta1[S1] <- 5
#' y1 <- X[, S1] %*% beta1[S1] + rnorm(n)

#' S2 <- 3:4 # second support
#' beta2 <- numeric(p)
#' beta2[S2] <- 5
#' y2 <- X[, S2] %*% beta2[S2] + rnorm(n)
#' 
#' 
#' ## Monte Carlo simulation giving empirical distribution for the statistic P (see paper below):
#' qut.MC.output <-  qut.MC(X, parallel = FALSE, MCrep = 50)
#' 
#' ## lasso-zero estimates:
#' 
#' ## for y1 with alpha = 0.1:
#' lass01 <- lass0(X, y1, alpha = 0.1, qut.MC.output = qut.MC.output, parallel = FALSE)
#' plot(lass01)
#' 
#' ## for y2 with alpha = 0.05:
#' lass02 <- lass0(X, y2, alpha = 0.05, qut.MC.output = qut.MC.output, parallel = FALSE)
#' plot(lass02)
#' 
#' @seealso  \code{\link{lass0}}
#'
#' @references Descloux, P., & Sardy, S. (2018). Model selection with
#'   lasso-zero: adding straw to the haystack to better find needles. arXiv
#'   preprint arXiv:1805.05133. \url{https://arxiv.org/abs/1805.05133}
#'
#' @references  Giacobino, C., Sardy, S., Diaz-Rodriguez, J., & Hengartner, N.
#'   (2017). Quantile universal threshold. Electronic Journal of Statistics,
#'   11(2), 4701-4722.
#'   

qut.MC <- function(X, q = nrow(X), M = 30, alpha = NULL, sigma = NULL,
                   intercept = TRUE, standardizeX = TRUE, standardizeG = NULL,
                   MCrep = 1.e2, GEVapprox = TRUE, parallel = FALSE) {
    
    ## check for some errors / warnings in the arguments
    if (M < 1) {
        stop("M must be a number >= 1")
    } else if (M > 1 & q == 0) {
        warning("q = 0 --> no noise dictionary is used, therefore M is set to 1")
        M <- 1
    } 
    
    lass0settings <- list(q = q, M = M, sigma = sigma, intercept = intercept,
                          standardizeX = standardizeX, standardizeG = standardizeG)
    
    ##
    X <- as.matrix(X)
    n <- nrow(X)
    if (intercept) X <- t(t(X) - colMeans(X))
    if (standardizeX) X <- t(t(X) / apply(X, 2, sd))
    if (is.null(standardizeG)) standardizeG <- standardizeX
    
    ## Monte Carlo simulation:
    if (parallel) {
        print("foreach loop starts")
        allMC <- foreach::foreach(r = 1:MCrep, .combine = "c", .packages = "lpSolve") %dopar% {
            print(paste0("MC rep", r))
            eps <- rnorm(n)
            if (intercept) eps <- eps - mean(eps)
            extBP <- MextBP(X = X, y = eps, q = q, M = M,
                            meancenterG = intercept, standardizeG = standardizeG,
                            parallel = FALSE, returnGammas = is.null(sigma))
            if (is.null(sigma)){
                max(abs(extBP$medBeta)) / extBP$madGammas
            }else{
                sigma * max(abs(extBP$medBeta))
            }
        }
    } else {
        allMC <- numeric(MCrep)
        for (r in 1:MCrep) {
            print(paste0("MC rep", r))
            eps <- rnorm(n)
            if (intercept) eps <- eps - mean(eps)
            extBP <- MextBP(X = X, y = eps, q = q, M = M,
                            meancenterG = intercept, standardizeG = standardizeG,
                            parallel = FALSE, returnGammas = is.null(sigma))
            if (is.null(sigma)) {
                allMC[r] <- max(abs(extBP$medBeta)) / extBP$madGammas
            } else {
                allMC[r] <- sigma * max(abs(extBP$medBeta))
            }
        }
    }
    
    ## if approximation by GEV distribution:
    if (GEVapprox) {
        print("fitting GEV distribution")
        GEVfit <- ismev::gev.fit(allMC)
        GEVpar <- GEVfit$mle
        if (!is.null(alpha)) {
            if (GEVpar[3] == 0) {
                upperQuant <- GEVpar[1] - GEVpar[2] * log(-log(1-alpha))
            } else {
                upperQuant <- GEVpar[1] + GEVpar[2] / GEVpar[3] * (1/(-log(1-alpha))^GEVpar[3] - 1)
            }
        } else {
            upperQuant <- NULL
        }
    } else {
        GEVpar <- "GEV parameters were not estimated"
        if (!is.null(alpha)) {
            upperQuant <- quantile(allMC, 1-alpha)
        } else {
            upperQuant <- NULL
        }
    }
    
    ## output:
    out <- list()
    out$allMC <- allMC
    out$GEVpar <- GEVpar
    out$upperQuant <- upperQuant
    out$call <- match.call()
    out$lass0settings <- lass0settings 
    class(out) <- "qut.MC"
    out
}
