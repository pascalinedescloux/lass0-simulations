generateX <- function (type = c("smallGauss", "wideGauss", "TV300", "riboflavin")) {
    
    # Generate the required design matrix X:
    # 
    # 
    # smallGauss: 100x200 i.i.d. Gaussian, mean-centered and standardized
    # wideGauss: 100x1000 i.i.d. Gaussian, mean-centered and standardized
    # TV300: matrix for segmentation problem with signal length of 300,
    #   (mean-centered, but not standardized)
    # riboflavin: design matrix from riboflavin dataset, mean-centered and standardized
    
    set.seed(2019)
    
    if (type == "smallGauss") {
        n <- 100
        p <- 200
        X <- matrix(rnorm(n*p), n, p)
    } else if (type == "wideGauss") {
        n <- 100
        p <- 1000
        X <- matrix(rnorm(n*p), n, p)
    } else if (type == "TV300") {
        n <- 300
        X <- matrix(0, n, n-1)
        X[col(X) < row(X)] <- 1
    } else if (type == "riboflavin") {
        library(qut)
        data(riboflavin)
        X <- riboflavin$x
        class(X) <- class(X)[-match("AsIs", class(X))]
    }
    
    X <- t(t(X) - colMeans(X))
    if (type != "TV300") {
        X <- t(t(X) / apply(X, 2, sd))
    }
    
    return(X)
}