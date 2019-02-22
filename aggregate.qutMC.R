## file to aggregated qut-MC results 
## (when performed in different joy arrays on baobab)

type <- "riboflavin"
arrays <- 1:10

for (k in arrays) {
    load(paste0("qutMC.", type, ".array", k, ".Rda"))
}

allLists <- mget(paste0("qutMC.", type, ".array", arrays))

newList <- list()
newList$allMC <- as.vector(sapply(allLists, function(L) L$allMC))
newList$GEVpar <-"GEV parameters were not estimated"
newList$lass0settings <- get(paste0("qutMC.", type, ".array", arrays[1]))$lass0settings
newList$arrays <- arrays
class(newList) <- "qut.MC"

assign(paste0("qutMC.", type), newList)
save(list=paste0("qutMC.", type), file=paste0("qutMC.", type, ".Rda"))


