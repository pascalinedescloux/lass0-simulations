ncv.byGIC <- function(X, y, penalty = c("SCAD", "MCP"), sigma, standardize = FALSE){
    library(ncvreg)
    n <- nrow(X)
    p <- ncol(X)
    GIC <- function(coefs){
        sum((y - coefs[1] - X %*% coefs[-1])^2) + sigma^2 * log(log(n)) * log(p) * sum(coefs[-1] != 0)
    }
    
    if(standardize){
        sdsX <- apply(X, 2, sd)
        X <- t(t(X) / sdsX)
    }else{sdsX <- rep(1, p)}
    
    res <- ncvreg(X, y, family = "gaussian", penalty = penalty)
    lambda.seq <- res$lambda
    all.GIC <- apply(res$beta, 2, GIC)
    ix <- which.min(all.GIC)
    lambda <- lambda.seq[ix]
    betahat <- res$beta[-1, ix]
    if(standardize) betahat <- betahat / sdsX
    betahat
}