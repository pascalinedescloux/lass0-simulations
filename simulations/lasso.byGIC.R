lasso.byGIC <- function(X, y, sigma, intercept=FALSE, standardize=FALSE){
    library(glmnet)
    n <- nrow(X)
    p <- ncol(X)
    if(intercept){
        meansX <- colMeans(X)
        X <- t(t(X) - meansX)
        meany <- mean(y)
        y <- y - meany
    }
    if(standardize){
        sdsX <- apply(X, 2, sd)
        X <- t(t(X) / sdsX)
    }else{sdsX <- rep(1, p)}
    
    res <- glmnet(X, y, standardize=FALSE, intercept=FALSE)
    lambda.seq <- res$lambda
    GIC <- sapply(1:length(lambda.seq), function(i){
        sum((y - X %*% res$beta[, i])^2) + sigma^2 * log(log(n)) * log(p) * sum(res$beta[, i] != 0)
    })
    ix <- which.min(GIC)
    betahat <- res$beta[, ix]
    betahat <- betahat / sdsX
    if(intercept) intercept <- meany - betahat %*% meansX
    
    out <- list()
    out$betahat <- betahat
    if(intercept) out$intercept <- intercept
    return(out)
}