### plots (with ggplot2):

library(ggplot2)
library(reshape2)

## load data:
type <- "TV300"
load(paste0(type, ".Rda"))
results <- get(type)
rm(list = type)

## to compute FWER, replace Nfalse by 1_{Nfalse > 0}
results[, dimnames(results)$criterion == "Nfalse", , ] <- as.numeric(results[, dimnames(results)$criterion == "Nfalse", , ] > 0)
dimnames(results)$criterion[dimnames(results)$criterion == "Nfalse"] <- "FWER"

## compute mean values
results <- apply(results, 1:3, mean)
s.values <- as.numeric(dimnames(results)$sparsity)
resnames <- dimnames(results)

## removing MX-knockoffs with offset = 1:
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

## long format:
Lresults <- melt(results, value.name = "mean")
Lresults$criterion <- gsub("tpp", "TPR", Lresults$criterion)
Lresults$criterion <- gsub("fdp", "FDR", Lresults$criterion)
Lresults$criterion <- factor(Lresults$criterion, levels = c("sign.recovery",
                                                            "supp.recovery",
                                                            "FWER",
                                                            "FDR", "TPR",
                                                            "modelsize"))
## colors and point styles:
colour.val <- c("lasso-zero" = "black",
                "lasso" = "red",
                "SCAD" = "olivedrab",
                "stability selection" = "orange",
                "knockoffs" = "blue",
                "SLOPE" = "magenta",
                "wild binary segmentation" = "purple")
colour.val <- colour.val[names(colour.val) %in% levels(Lresults$estimator)]
shape.val <- c("lasso-zero" = 1,
               "lasso" = 2,
               "SCAD" = 3,
               "stability selection" = 4,
               "knockoffs" = 5,
               "SLOPE" = 6,
               "wild binary segmentation" = 0)
# shape.val <- c("lasso-zero" = 19,
#                 "lasso" = 24,
#                 "SCAD" = 13,
#                 "stability selection" = 7,
#                 "knockoffs" = 23,
#                 "SLOPE" = 25,
#                 "wild binary segmentation" = 15)
shape.val <- shape.val[names(shape.val) %in% levels(Lresults$estimator)]

# support recovery:
ggplot(Lresults[Lresults$criterion == "supp.recovery", ], 
       aes(x = sparsity, y = mean)) +
    geom_point(aes(colour = estimator,  shape = estimator)) +
    geom_line(aes(colour = estimator), size = 0.2) +
    scale_colour_manual(values = colour.val) +
    scale_shape_manual(values = shape.val) +
    labs(colour = "", shape = "") +
    ylim(0, 1) +
    ylab(expression(P(hat(S) == S^0))) + 
    xlab(expression(s^0)) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal", 
          legend.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.background = element_rect(fill = "gray95")) +
    guides(col = guide_legend(nrow=2))



# TPR-FDR-FWER:
ggplot(Lresults[Lresults$criterion %in% c("TPR", "FDR", "FWER"), ],
       aes(x = sparsity, y = mean)) +
    facet_grid(criterion ~ .) +
    geom_point(aes(colour = estimator,  shape = estimator)) +
    geom_line(aes(colour = estimator), size = 0.2) +
    geom_hline(data = data.frame(criterion=c("FWER", "FDR"), y=rep(0.05, 2)), 
               mapping=aes(yintercept=y), linetype = "dashed") +
    scale_colour_manual(values = colour.val) +
    scale_shape_manual(values = shape.val) +
    labs(colour = "", shape = "") +
    ylim(0, 1) +
    ylab("") + xlab(expression(s^0)) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal", 
          legend.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.background = element_rect(fill = "gray95")) +
    guides(col = guide_legend(nrow=2))

