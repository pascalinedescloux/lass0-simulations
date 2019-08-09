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
