
### from nsRFA
### need change for the precision of the P value 
### also uses lmom functions rather than the nsRFA ones now

#' Modified functions from nsRFA
#'
#' @param x the vector to be tested
#' @param Nsim the number of simulations which should be run
#' @return The A2 statistics and the CORRECT p value
#' @aliases gofGLOtest
#' @export
gofGEVtest_ip <- function (x, Nsim = 10000) 
{
    x <- sort(x)
    n <- length(x)
    Lmom.x <- samlmu(x)
    par <- pelgev(c(Lmom.x["l_1"], Lmom.x["l_2"], Lmom.x["t_3"]))
    F <- suppressWarnings(cdfgev(x,c(par["xi"], par["alpha"], par["k"])))
    F[F < 1e-08] <- 1e-08
    F[F > 0.99999999] <- 0.99999999
    A2 <- -n - (1/n) * sum((seq(1, 2 * n - 1, by = 2)) * log(F) + 
        (seq(2 * n - 1, 1, by = -2)) * log(1 - F))
    A2s <- rep(NA, Nsim)
    for (i in 1:Nsim) {
        x.sim <- quagev(runif(n),c(par["xi"], par["alpha"], par["k"]))
        x.sim <- sort(x.sim)
        Lmom.xsim <- samlmu(x.sim)
        par.sim <- pelgev(c(Lmom.xsim["l_1"], Lmom.xsim["l_2"], Lmom.xsim["t_3"]))
        F <- suppressWarnings(cdfgev(x.sim,c(par.sim["xi"], par.sim["alpha"], par.sim["k"])))
        F[F < 1e-08] <- 1e-08
        F[F > 0.99999999] <- 0.99999999
        A2s[i] <- -n - (1/n) * sum((seq(1, 2 * n - 1, by = 2)) * 
            log(F) + (seq(2 * n - 1, 1, by = -2)) * log(1 - F))
    }
    ecdfA2s <- ecdf(A2s)
    probabilita <- 1 - ecdfA2s(A2)
    output <- c(A2, probabilita)
    names(output) <- c("A2", "P")
    return(output)
}


gofGLOtest_ip <- function (x, Nsim = 10000) {
    x <- sort(x)
    n <- length(x)
    Lmom.x <- samlmu(x)
    par <- pelglo(c(Lmom.x["l_1"], Lmom.x["l_2"], Lmom.x["t_3"]))
    F <- suppressWarnings(cdfglo(x,c(par["xi"], par["alpha"], par["k"])))
    F[F < 1e-08] <- 1e-08
    F[F > 0.99999999] <- 0.99999999
    A2 <- -n - (1/n) * sum((seq(1, 2 * n - 1, by = 2)) * log(F) + 
        (seq(2 * n - 1, 1, by = -2)) * log(1 - F))
    A2s <- rep(NA, Nsim)
    for (i in 1:Nsim) {
        x.sim <- quaglo(runif(n),c(par["xi"], par["alpha"], par["k"]))
        x.sim <- sort(x.sim)
        Lmom.xsim <- samlmu(x.sim)
        par.sim <- pelglo(c(Lmom.xsim["l_1"], Lmom.xsim["l_2"], Lmom.xsim["t_3"]))
        F <- suppressWarnings(cdfglo(x.sim,c(par.sim["xi"], par.sim["alpha"], par.sim["k"])))
        F[F < 1e-08] <- 1e-08
        F[F > 0.99999999] <- 0.99999999
        A2s[i] <- -n - (1/n) * sum((seq(1, 2 * n - 1, by = 2)) * 
            log(F) + (seq(2 * n - 1, 1, by = -2)) * log(1 - F))
    }
    ecdfA2s <- ecdf(A2s)
    probabilita <- 1 - ecdfA2s(A2)
    output <- c(A2, probabilita)
    names(output) <- c("A2", "P")
    return(output)
}

