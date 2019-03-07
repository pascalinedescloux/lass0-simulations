## aggregate result of simulation for different sparsity indices:

library(abind)

#######################
type <- "wideGauss"
s.values <- 0:15
#######################

for (s in s.values) {
    load(paste0(type, ".sparsity", s, ".Rda"))
}

newArray <- abind(mget(paste0(type, ".sparsity", s.values)), along = 0)
names(dimnames(newArray)) <- c("sparsity", "criterion", "estimator", "rep")
dimnames(newArray)$rep <- 1:(dim(newArray)[length(dim(newArray))])
dimnames(newArray)$sparsity <- s.values

assign(type, newArray)
save(list = type, file = paste0(type, ".Rda"))

