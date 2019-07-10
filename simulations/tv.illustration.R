## TV illustration

library(lass0)
library(doParallel)
registerDoParallel(cores=detectCores())
source("generateX.R")
source("comp.TV300.Gnorm.R")

# design
X <- generateX(type = "TV300")
n <- nrow(X)
p <- ncol(X)
# standardization for noise dictionaries:
N <- 1.e4
standardizeG <- comp.TV300.Gnorm(N = N, q = n, alpha = 0.05, 
                                 parallel = TRUE, seeds = 1:N) 
# sparsity:
s <- 3
# locations of nonzero coefficients
S <- (1:s) * round(n/(s+1))
# amplitude of nonzero coefficients:
amp <- 2.5
# noise level
sigma <- 1


set.seed(8)

# piecewise constant function (jumps have same sign)
#jumps <- runif(s, amp.min, amp.max)
jumps <- rep(amp, s)
f <- rep(cumsum(c(0, jumps)), times=diff(c(0,S, n)))
y <- f + rnorm(n, 0, sigma)
yc <- y - mean(y)

# BPsol: (q, M) = (0, 1)
BPsol <- lass0(X, yc, tau = 0, q = 0, M = 1, standardizeX = FALSE, 
               standardizeG = standardizeG, ols = FALSE)$coefficients
# with 1 noise dictionary (q, M) = (n, 1)
M1sol <- lass0(X, yc, tau = 0, q = n, M = 1, standardizeX = FALSE, 
               standardizeG = standardizeG, ols = FALSE)$coefficients
# with 30 noise dictionaries: (q, M) = (n, 30)
M30sol <- lass0(X, yc, tau = 0, q = n, M = 30, standardizeX = FALSE, 
               standardizeG = standardizeG, ols = FALSE)$coefficients



par(mar=c(5, 6, 2, 0.5) + 0.1)
#ylab <- expression(hat(beta)[j]^"\u2113"[1])
ylab <- expression(hat(beta)[j]^"med")

## BP:
plot(BPsol, col=1+(1:p)%in%S,
     xlab="j", ylab="", type="h", cex.lab=1.5)
segments(x0=S, y0=rep(0, s), x1=S, y1=BPsol[S],
         col=2, lwd=2)
mtext(ylab, side=2, line=2, las=2, cex=2, at=1)

## M = 1:
plot(M1sol, col=1+(1:p)%in%S,
     xlab="j", ylab="", type="h", cex.lab=1.5)
segments(x0=S, y0=rep(0, s), x1=S, y1=M1sol[S],
         col=2, lwd=2)
mtext(ylab, side=2, line=2, las=2, cex=2, at=1)

## M = 30:
plot(M30sol, col=1+(1:p)%in%S,
     xlab="j", ylab="", type="h", cex.lab=1.5)
segments(x0=S, y0=rep(0, s), x1=S, y1=M30sol[S],
         col=2, lwd=2)
mtext(ylab, side=2, line=2, las=2, cex=2, at=1)

par(mar=c(5, 4, 4, 2) + 0.1)
