## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.retina=2, fig.width=7, fig.height=5)

## ----plotSyst, fig.show='asis'-------------------------------------------
par(mar = c(2, 2, 2, 0.5))
syst <- data.frame(Year = seq(1976,2015), 
                   Flow = c(300.346, 527.381, 410.964, 224.615, 370.688, 275.477, 
                            345.358, 439.179, 117.119, 371.371, 386.094, 458.024, 
                            286.608, 264.741, 262.419, 453.842, 715.718, 480.023, 
                            316.092, 479.976, 457.083, 317.254, 378.389, 249.409,
                            384.964, 228.965, 291.549, 370.726, 235.923, 229.726, 
                            376.540, 326.176, 553.550, 318.610, 105.712, 309.667, 
                            375.903, 523.235, 257.639, 229.343))
with(syst, plot(Year,Flow,type="h"))

## ------------------------------------------------------------------------
lmom::samlmu(syst$Flow) ## at-site L-moments
## GLO parameter estimates via L-moment approach
systOnlyLmom <- lmom::pelglo(lmom::samlmu(syst$Flow)) 
## GLO parameter estimates via Maximum Likelihood approach
systOnlyMaxL <- ilaprosUtils::glo.fit(syst$Flow, show = FALSE)
systOnlyMaxL$conv ## if 0 model has converged
systOnlyLmom; systOnlyMaxL$mle ## not very different

## ----plotRetSyst, fig.show='asis'----------------------------------------
## ?ilaprosUtils::retPlot
rlSOML <- ilaprosUtils::retPlot(systOnlyMaxL, pch = 16, sign.alpha = 0)
lines(log(rlSOML$p/(1-rlSOML$p)), 
      lmom::quaglo(f = rlSOML$p, para = systOnlyLmom), col = 2)

## ----plothist, fig.show='asis'-------------------------------------------
hist <- data.frame(Year = c(1770, 1850, 1920),  ## only needed for plotting
                   Flow = c(840, 820, 870))
with(syst, plot(Year, Flow, type="h", xlim =  c(1976-250, 2015), ylim = c(0,900)))
lines(hist$Year, hist$Flow, type="h", col = 2)
abline(h = 800, col = 4, lty = 4) ## perception threshold

## ------------------------------------------------------------------------
ll <- ilaprosUtils::lmom.hist.fit(
            xdat = c(hist$Flow,syst$Flow), ## all data combined
            h = 250, ## historical period length
            X0 = 800, ## perception threshold
            k = 3, ## number of historical peaks
            nmom = 4) ## number of L-mom parameters 
        
## GLO parameter estimates via L-moment approach
withHistLmom <- lmom::pelglo(ll); rm(ll)
## GLO parameter estimates via Maximum Likelihood approach
withHistMaxL <- ilaprosUtils::glo.hist.fit(
            xdat = c(hist$Flow,syst$Flow), ## all data combined
            h = 250, ## historical period length
            X0 = 800, ## perception threshold
            k = 3, ## number of historical peaks
            show = FALSE)
withHistMaxL$conv ## if 0 model has converged - check!
withHistLmom; withHistMaxL$mle ## somewhat different

## ----retHist, fig.show='asis'--------------------------------------------
## ?ilaprosUtils::retPlot
## retPlot recognises that some data are from an historical record
## displayed as filled squares
rlWHML <- ilaprosUtils::retPlot(withHistMaxL, 
                p = c(seq(0.01,0.92,l=50),seq(0.92,.998,l=20)),
                pch = 16, sign.alpha = 0, lty = 2)
## Sysy. only  max lik
lines(log(rlWHML$p/(1-rlWHML$p)), 
      lmom::quaglo(f = rlWHML$p, para = systOnlyMaxL$mle), col = 1)
## Sysyonly - Lmoments
lines(log(rlWHML$p/(1-rlWHML$p)), 
      lmom::quaglo(f = rlWHML$p, para = systOnlyLmom), col = 2)
## With Historical  info - Lmoments
lines(log(rlWHML$p/(1-rlWHML$p)), 
      lmom::quaglo(f = rlWHML$p, para = withHistLmom), col = 2, lty = 2)

legend("topleft", col = c(1,1,2,2), lty = c(1,2,1,2), bty = "n",
       legend = c("Syst. only - ML", "With Hist. - ML",
                  "Syst. only - Lmom", "With Hist. - Lmom"))
legend("bottomright", pch = c(16,15), bty = "n",
       legend = c("Systematic data", "Historical data"))

## ------------------------------------------------------------------------
## GLO parameter estimates via Maximum Likelihood approach
## Under Binomial censoring
withBinCMaxL <- ilaprosUtils::glo.hist.fit(
            ### The first points only need to be larger than X0           
            xdat = c(801, 801, 801 ,syst$Flow), ## all data combined
            h = 250, ## historical period length
            X0 = 800, ## perception threshold
            k = 3, ## number of historical peaks
            binomialcens = TRUE,
            show = FALSE)
withBinCMaxL$conv ## if 0 model has converged - check!
withBinCMaxL$mle; withHistMaxL$mle ## little difference

