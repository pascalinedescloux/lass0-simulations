BP <- function(X, y) {
    p <- ncol(X)
    out.lp <- lpSolve::lp(objective.in = rep(1, 2*p), const.mat = cbind(X, -X), 
                          const.dir = "==", const.rhs =y )
    if (out.lp$status == 2) stop("linear system y = X beta + G gamma has no solution; choose larger q (e.g. q >= nrow(X))" )
    # TO IMPLEMENT: LS solution with minimal l1 norm ?
    sol <- out.lp$solution[1:p] - out.lp$solution[-(1:p)]
    sol
}

MextBP <- function(X, y, q = ncol(X), M = 30, meancenterG, standardizeG, 
                    parallel = FALSE, returnGammas = TRUE) {
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(q)) q <- n
    if (parallel) {
        BPsols <- foreach::foreach(m = 1:M, .combine = "cbind", .packages = "lpSolve") %dopar% {
            G <- matrix(rnorm(n*q), n, q)
            if (meancenterG) G <- t(t(G) - colMeans(G))
            if (is.numeric(standardizeG)) {
                G <- standardizeG * t(t(G) / apply(G, 2, function(v) sqrt(sum(v^2))))
            } else if (standardizeG) {
                G <- t(t(G) / apply(G, 2, sd))
            }
            BP(cbind(X, G), y)
        }
    } else {
        BPsols <- array(NA, dim=c(p+q, M))
        for (m in 1:M) {
            G <- matrix(rnorm(n*q), n, q)
            if (meancenterG) G <- t(t(G) - colMeans(G))
            if (is.numeric(standardizeG)) {
                G <- standardizeG * t(t(G) / apply(G, 2, function(v) sqrt(sum(v^2))))
            } else if (standardizeG) {
                G <- t(t(G) / apply(G, 2, sd))
            }
            BPsols[, m] <- BP(cbind(X, G), y)
        }
    }
    
    Betas <- BPsols[1:p, ]
    if (M == 1) {
        medBeta <- Betas
    } else {
        medBeta <- apply(Betas, 1, median)
    }
    if (returnGammas) {
        Gammas <- BPsols[-(1:p), ]
        gammas <- as.vector(Gammas)
        madGammas <- mad(gammas[gammas != 0])
    } else {
        Gammas <- NULL
        madGammas <- NULL
    }
    list(Betas = Betas, medBeta = medBeta, 
         Gammas = Gammas, madGammas = madGammas)
}