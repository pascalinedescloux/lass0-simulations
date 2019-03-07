## lasso path


library(glmnet)
library(lass0)
library(MASS)
source("generateX.R")

set.seed(1)

# setting: wideGauss(100x1000)
X <- generateX(type = "wideGauss")
n <- nrow(X)
p <- ncol(X)
s <- 10
amp <- 0.75
sigma <- 1

# signal
beta <- rep(0, p)
beta[1:s] <- sign(rnorm(s)) * amp
y <- X %*% beta + rnorm(n, 0, sigma)
y <- y - mean(y)

# IR:
S <- 1:s
signs <- sign(beta[S])
v <- t(X[, -S]) %*% X[, S] %*% solve(t(X[, S]) %*% X[, S]) %*% signs
as.numeric(max(abs(v)) < 1)


# BPsol:
BPsol <- lass0(X, y, tau = 0, q = 0, M = 1, intercept = FALSE, 
               standardizeX = FALSE, ols = FALSE)$coefficients
M <- max(abs(BPsol[-(1:s)]))

# glmnet
res <- glmnet(X, y, standardize=FALSE, intercept=FALSE, lambda.min.ratio = 0.0001)

# plot
col.ix <- rep(1, p)
col.ix[1:s] <- 2

par(mar=c(4, 8, 1, 4) + 0.1)
xlab <- expression(lambda)
ylab <- expression(hat(beta)[j]^lasso)
linetype <- rep(1:6, p)
linetype[S] <- rep(1, s)
matplot(c(res$lambda, 0), rbind(t(res$beta), BPsol), type="l", col=col.ix, lty=linetype,
        xlab=xlab, ylab="", cex.lab=2, cex.axis=1.5, cex.main=2)
mtext(ylab, 2, line=3, las=1, cex=2)
mtext(expression(tau), 4, line=2, at=1.1*M, las=1, cex=1.5)
mtext(expression(-tau), 4, line=2, at=-1.1*M, las=1, cex=1.5)
abline(h=c(1.1 * M, -1.1*M), lty=2, col=4)
par(mar=c(5, 4, 4, 2) + 0.1)

