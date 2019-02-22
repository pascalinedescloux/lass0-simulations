## plots:

type <- "riboflavin"
load(paste0(type, ".Rda"))
results <- get(type)
rm(list=type)
results <- apply(results, 1:3, mean)
s.values <- as.numeric(dimnames(results)$sparsity)
resnames <- dimnames(results)
# removing MX-knockoffs with offset = 1:
results <- results[ , , resnames$estimator != "MXknockoffs.off1"]
# modify estimator names:

resnames <- dimnames(results)
pch.ix <- 1:length(resnames$estimator)
col.ix <- 1:length(resnames$estimator)
xlab <- expression(s^0)

# support recovery:
par(mar=c(5, 6, 1, 1) + 0.1)
matplot(s.values, results[, resnames$criterion == "supp.recovery", ],
        type = "b", pch = pch.ix, lty = 1, main = "", xlab = xlab, 
        ylab = expression(P(hat(S) == S^0)), col = col.ix,
        cex.lab=2, cex.axis=1.5, cex.main=2)
if (type == "smallGauss") {
    legend(x = 12, y = 1, legend = resnames$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
} else if (type == "wideGauss") {
    legend(x = 10, y = 1, legend = resnames$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
}
par(mar=c(5, 4, 4, 2) + 0.1)

# fdr & tpr
par(mfrow=c(2, 1), mar=c(5, 5, 1, 2) + 0.1)
matplot(s.values, results[ , resnames$criterion == "fdp", ], type = "b", 
        pch=pch.ix, lty=1, xlab="", ylab="", col=col.ix, cex.lab=2, cex.axis=1.5, cex.main=2)
title(ylab = "FDR", line = 3, cex.lab = 2)
abline(h=0.05, lty=2)
par(mar=c(c(5, 5, 1, 2) + 0.1))
matplot(s.values, results[ , resnames$criterion == "tpp", ], type="b", 
        pch=pch.ix, lty=1, xlab=xlab, ylab="", col=col.ix, cex.lab=2, cex.axis=1.5, cex.main=2)
title(ylab="TPR", line=3, cex.lab=2)
if (type == "smallGauss") {
    legend(x = 0, y = 0.6, legend = resnames$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.75) 
} else if (type == "wideGauss") {
    legend(x = 0, y = 0.6, legend = resnames$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.75) 
}
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
