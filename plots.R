## plots:

type <- "TV300"
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
est.names <- resnames$estimator
est.names <- gsub("lassoZero", "lasso-zero", est.names)
est.names <- gsub("MXknockoffs.off0", "knockoffs", est.names)
est.names <- gsub("SLOPE.gauss", "SLOPE", est.names)
est.names <- gsub("stabilitySelection", "stability selection", est.names)
est.names <- gsub("wildBinarySeg", "wild binary segmentation", est.names)
dimnames(results)$estimator <- est.names

# colors and point styles:
all.est <- c("lasso-zero", "lasso", "SCAD", "stability selection", 
             "knockoffs", "SLOPE", "wild binary segmentation")
all.col <-  c("black", "red", "limegreen", "orange", "blue", "magenta", "purple")
pch.ix <- match(dimnames(results)$estimator, all.est)
col.ix <- all.col[match(dimnames(results)$estimator, all.est)]

# x-label
xlab <- expression(s^0)

# support recovery:
par(mfrow=c(1,1), mar=c(5, 6, 1, 1) + 0.1)
matplot(s.values, results[, resnames$criterion == "supp.recovery", ],
        type = "b", pch = pch.ix, lty = 1, main = "", xlab = xlab, 
        ylab = expression(P(hat(S) == S^0)), col = col.ix,
        cex.lab=2, cex.axis=1.5, cex.main=2)
if (type == "smallGauss") {
    legend(x = 15, y = 1, legend = resnames$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
} else if (type == "wideGauss") {
    legend(x = 8, y = 1, legend = resnames$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
} else if (type == "riboflavin") {
    legend(x = 6, y = 1, legend = resnames$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
} else if(type == "TV300") {
    legend(x=5, y = 1, legend = resnames$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
}
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

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
} else if (type == "riboflavin") {
    legend(x = 0, y = 0.5, legend = resnames$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.75) 
} else if (type == "wideGauss") {
    legend(x = 0, y = 0.6, legend = resnames$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.75) 
} else if (type == "TV300") {
    legend(x = 0, y = 0.8, legend = resnames$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.75) 
}
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
