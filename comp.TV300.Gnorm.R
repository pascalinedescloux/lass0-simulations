comp.TV300.Gnorm <- function(N = 1.e4, q, alpha, parallel = TRUE,
                             seeds = NULL) {
    
    # For the TV300 design matrix X, computes the constant C = TV300.Gnorm
    # such that supnorm(X^t eps) and supnorm(Gnorm^t eps) have the same upper alpha
    # quantile, where Gnorm is a n x q Gaussian matrix that has been mean-centered
    # and standardized so that its euclidean norm is C*G.
    
    # N: sample size for the empirical distribution
    # q: size of noise dictionaries. By default: q = nrow(X)
    # alpha: upper quantile at which we want to align both statistics for X and G
    # parallel: if TRUE, realizations are generated in parallel with foreach
    # seeds: (for reproducibility) if parallel = TRUE, vectors of N different seeds,
    #       if parallel =  FALSE, a numerical value.
    
    source("generateX.R")
    X <- generateX(type = "TV300")
    n <- nrow(X)
    if (is.null(q)) q <- n
    
    if (parallel) {
        library(doParallel)
        registerDoParallel(cores = detectCores())
        supnorms <- foreach(m = 1:N, .combine = "cbind") %dopar% {
            print(m)
            if(!is.null(seeds)) set.seed(seeds[m])
            eps <- rnorm(n)
            G <- matrix(rnorm(n*q), n, q)
            G <- t(t(G) - colMeans(G))
            G <- t(t(G) / apply(G^2, 2, function(g) sqrt(sum(g))))
            c(max(abs(t(G) %*% eps)), max(abs(t(X) %*% eps)))
        }
        TV300.Gnorm <- quantile(supnorms[2, ], 1-alpha) / quantile(supnorms[1, ], 1-alpha)
    } else { 
        if(!is.null(seeds)) set.seed(seeds)
        comp.stat <- function(eps) {
            G <- matrix(rnorm(n*q), n, q)
            G <- t(t(G) - colMeans(G))
            G <- t(t(G) / apply(G^2, 2, function(g) sqrt(sum(g))))
            max(abs(t(G) %*% eps))
        }
        Eps <- matrix(rnorm(N*n), n, N)
        supnormsG <- apply(Eps, 2, comp.stat)
        supnormsX <- apply(abs(t(X) %*% Eps), 2, max)
        TV300.Gnorm <- quantile(supnormsX, 1-alpha) / quantile(supnormsG, 1-alpha)
    }
    
    
    return(TV300.Gnorm)
}