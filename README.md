# Lasso-Zero
## Pascaline Descloux, Sylvain Sardy

R package and simulations for Lasso-Zero, a new variable selection method for (high-dimensional) linear regression.
See https://arxiv.org/abs/1805.05133.

### R package:
To install the R package lass0 directly from GitHub, run the following in R:

```{r}
library(devtools)
install_github(repo="pascalinedescloux/lasso-zero", subdir="lass0")
```

### Simulations:
Contains the code that generated all figures of the paper (adapted for running on the Baobab cluster of University of Geneva). 
The main script simuSparsityFixed.R runs a simulation for fixed design matrix and sparsity index s.
