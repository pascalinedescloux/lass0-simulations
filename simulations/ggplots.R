## plots (with ggplot2):

library(ggplot2)
library(reshape2)

type <- "smallGauss"
load(paste0(type, ".Rda"))
results <- get(type)
rm(list = type)
results <- apply(results, 1:3, mean)
s.values <- as.numeric(dimnames(results)$sparsity)
resnames <- dimnames(results)
# removing MX-knockoffs with offset = 1:
results <- results[ , , resnames$estimator != "MXknockoffs.off1"]

# modify estimator names:
est.names <- dimnames(results)$estimator
est.names <- gsub("lassoZero", "lasso-zero", est.names)
est.names <- gsub("MXknockoffs.off0", "knockoffs", est.names)
est.names <- gsub("knockoff.Xfixed.off0", "knockoffs", est.names)
est.names <- gsub("SLOPE.gauss", "SLOPE", est.names)
est.names <- gsub("stabilitySelection", "stability selection", est.names)
est.names <- gsub("wildBinarySeg", "wild binary segmentation", est.names)
dimnames(results)$estimator <- est.names

# long format:
Lresults <- melt(results, value.name = "mean")
Lresults$criterion <- gsub("tpp", "TPR", Lresults$criterion)
Lresults$criterion <- gsub("fdp", "FDR", Lresults$criterion)
Lresults$criterion <- factor(Lresults$criterion, levels = c("sign.recovery",
                                                            "supp.recovery",
                                                            "FDR", "TPR",
                                                            "modelsize",
                                                            "Nfalse"))
# colors and point styles:
all.est <- c("lasso-zero", "lasso", "SCAD", "stability selection", 
             "knockoffs", "SLOPE", "wild binary segmentation")
all.col <-  c("black", "red", "limegreen", "orange", "blue", "magenta", "purple")
pch.ix <- match(dimnames(results)$estimator, all.est)
col.ix <- all.col[match(dimnames(results)$estimator, all.est)]

# support recovery:
ggplot(Lresults[Lresults$criterion == "supp.recovery", ], 
       aes(x=sparsity, y=mean, colour = estimator, shape = estimator)) +
    geom_line() +
    geom_point() +
    scale_colour_manual(values=c("lasso-zero" = "black",
                                 "lasso" = "red",
                                 "SCAD" = "limegreen",
                                 "stability selection" = "orange",
                                 "knockoffs" = "blue",
                                 "SLOPE" = "magenta",
                                 "wild binary segmentation" = "purple")) +
    scale_shape_manual(values = c("lasso-zero" = 1,
                                  "lasso" = 2,
                                  "SCAD" = 3,
                                  "stability selection" = 4,
                                  "knockoffs" = 5,
                                  "SLOPE" = 6,
                                  "wild binary segmentation" = 0)) +
    labs(colour = "", shape = "") +
    ylab("") + xlab(expression(s^0)) +
    theme(legend.position = "bottom",
          legend.direction = "vertical", 
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 16))

# TPR-FDR:
ggplot(Lresults[Lresults$criterion %in% c("TPR", "FDR"), ], 
       aes(x = sparsity, y = mean, colour = estimator, shape = estimator)) +
    geom_line() +
    geom_point() +
    facet_grid(criterion ~ .) +
    scale_colour_manual(values=c("lasso-zero" = "black",
                                 "lasso" = "red",
                                 "SCAD" = "limegreen",
                                 "stability selection" = "orange",
                                 "knockoffs" = "blue",
                                 "SLOPE" = "magenta",
                                 "wild binary segmentation" = "purple")) +
    scale_shape_manual(values = c("lasso-zero" = 1,
                                  "lasso" = 2,
                                  "SCAD" = 3,
                                  "stability selection" = 4,
                                  "knockoffs" = 5,
                                  "SLOPE" = 6,
                                  "wild binary segmentation" = 0)) +
    labs(colour = "", shape = "") +
    ylab("") + xlab(expression(s^0)) +
    theme(legend.position = "bottom",
          legend.direction = "vertical", 
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 16))








# x-label
xlab <- expression(s^0)

# support recovery:
par(mfrow=c(1,1), mar=c(4, 6, 4, 1) + 0.1)
matplot(s.values, results[, resnames$criterion == "supp.recovery", ],
        type = "b", pch = pch.ix, lty = 1, main = "", xlab = xlab, 
        ylim = c(0, 1),
        ylab = expression(P(hat(S) == S^0)), col = col.ix,
        cex.lab=2, cex.axis=1.5, cex.main=2)
if (type == "smallGauss") {
    legend(x = 15, y = 1, legend = dimnames(results)$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
} else if (type == "wideGauss") {
    legend(x = 8, y = 1, legend = dimnames(results)$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
} else if (type == "riboflavin") {
    legend(x = 6, y = 1, legend = dimnames(results)$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
} else if(type == "TV300") {
    legend(x=5, y = 1, legend = dimnames(results)$estimator, pch = pch.ix,
           lty = 1, col = col.ix, cex=1.3)
}
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

# fdr & tpr
par(mfrow=c(2, 1), mar=c(4, 5, 1, 2) + 0.1)
matplot(s.values, results[ , resnames$criterion == "fdp", ], type = "b", 
        ylim = c(0, max(results[ , resnames$criterion == "fdp", ])),
        pch=pch.ix, lty=1, xlab="", ylab="", col=col.ix, cex.lab=2, cex.axis=1.5, cex.main=2)
title(ylab = "FDR", line = 3, cex.lab = 2)
abline(h=0.05, lty=2)
par(mar=c(c(4, 5, 1, 2) + 0.1))
matplot(s.values, results[ , resnames$criterion == "tpp", ], type="b", 
        ylim = c(0, max(results[ , resnames$criterion == "tpp", ], na.rm = TRUE)),
        pch=pch.ix, lty=1, xlab=xlab, ylab="", col=col.ix, cex.lab=2, cex.axis=1.5, cex.main=2)
title(ylab="TPR", line=3, cex.lab=2)
if (type == "smallGauss") {
    legend(x = 0, y = 0.6, legend = dimnames(results)$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.75) 
} else if (type == "riboflavin") {
    legend(x = 0, y = 0.5, legend = dimnames(results)$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.75) 
} else if (type == "wideGauss") {
    legend(x = 0, y = 0.6, legend = dimnames(results)$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.75) 
} else if (type == "TV300") {
    legend(x = 7, y = 0.5, legend = dimnames(results)$estimator, col = col.ix, 
           lty = 1, pch = pch.ix, cex = 0.5) 
}
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)