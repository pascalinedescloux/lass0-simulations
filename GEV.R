###  Approximating the upper alpha quantile of P using a GEV fit:

library(ismev) # for gev.fit function

## simulation setting:
type <- "riboflavin" # "smallGauss", "wideGauss", "riboflavin", "TV300"

## alpha value:
alpha <- 0.01

## loading results of the (large scale - M = 10,000) Monte Carlo simulation
load(paste0("qutMC.", type, ".Rda"))
x <- get((paste0("qutMC.", type)))$allMC
Pquantile <- quantile(x, probs=1-alpha)
# removing outliers 
# (in our results: exactly one for smallGauss and wideGauss, 
# none for riboflavin and TV300):
x <- x[x < 1.e9]

## GEV quantile function
Q <- function(p, mu, sigma, xi) {
    if (xi == 0) {
        mu - sigma * log(-log(p))
    } else {
        mu + sigma / xi * (1/(-log(p))^xi - 1)
    }
}

## subsamples:
R <- 200 # number of subsamples
m <- 100 # subsample size

## repeat quantile estimation with GEV
set.seed(1)
GEVquantiles <- numeric(R)
MCquantiles <- numeric(R)
for(r in 1:R){
    trainix <- sample(1:length(x), m)
    xtrain <- x[trainix]
    ## m-MC quantile:
    MCquantiles[r] <- quantile(xtrain, probs = 1-alpha)
    ## GEV fit:
    GEVfit <- gev.fit(xtrain)
    GEVpars <- GEVfit$mle
    ## quantiles:
    GEVquantiles[r] <- Q(1-alpha, mu = GEVpars[1], sigma = GEVpars[2], xi = GEVpars[3])
}

# boxplot(t(estQuantiles))
# abline(h=Pquantile)
# boxplot(t(abs(estQuantiles-Pquantile)/Pquantile))


# ## qqplot:
# GEVfit <- gev.fit(x)
# GEVpars <- GEVfit$mle
# GUMfit <- gum.fit(x)
# GUMpars <- GUMfit$mle

# prob <- seq(1.e-4, 1-1.e-4, length=length(x))
# plot(quantile(x, prob), Q(prob, mu=GEVpars[1], sigma=GEVpars[2], xi=GEVpars[3]), main="GEV QQplot")
# abline(a=0, b=1, lty=2)
# plot(quantile(x, prob), Q(prob, mu=GUMpars[1], sigma=GUMpars[2], xi=0), main="Gumbel QQplot")
# abline(a=0, b=1, lty=2)


## plots:
par(mfrow=c(1, 2))
par(mar=c(5, 4, 4, 5))
boxplot(cbind(GEVquantiles, MCquantiles), main=expression(bold("estimates of" ~ q[alpha])),
        names=c("", ""),
        cex.lab=2, cex.axis=1.5, cex.main=2)
mtext(text=c(expression(hat(q)[alpha]^GEV), expression(hat(q)[alpha]^MC[100])),
      side=1, line=3, at=1:2, cex=2)
abline(h=Pquantile, lty=2)
mtext(text=expression(hat(q)[alpha]^MC[10000]), side=4, line=0.5, cex=2,
      las=1, at=Pquantile)
prob <- seq(1.e-4, 1-1.e-4, length=length(x))
# fit GEV on all 10'000:
par(mar=c(5, 7, 4, 2))
GEVfit <- gev.fit(x)
GEVpars <- GEVfit$mle
plot(quantile(x, prob), Q(prob, mu=GEVpars[1], sigma=GEVpars[2], xi=GEVpars[3]), 
     main="Q-Q plot", xlab="empirical quantiles", ylab="GEV quantiles",
     cex.lab=2, cex.axis=1.5, cex.main=2)
abline(a=0, b=1, lty=2)
par(mfrow=c(1,1))
