test_that("distribution functions work", {
  expect_warning(object = dgev(c(1,2,4),3,c(2,3),0.2), 3)
  expect_length(object = dgev(c(1,2,4),3,3,0.2), 3)
  expect_length(object = pgev(2,3,c(2,3),0.2), 2)
  expect_length(object = rgev(n = 10,6,c(0.01,3),0.2), 10)
  expect_length(object = qgev(c(.2,.3,.6),6,c(1,2,3,1,2,3),0.2), 6)
  expect_length(object = dglo(c(1,2,4),3,3,0.2), 3)
  expect_length(object = pglo(c(2,3),2,3,0.2), 2)
  expect_length(object = rglo(n = 10,c(6,60),03,0.2), 10)
  expect_length(object = qglo(c(.2,.3,.6),6,c(1,2,3,1,2,3),0.2), 6)
  expect_length(object = dkappa(4,6,c(1,2,3),0.2,-0.3), 3)
  expect_length(object = pkappa(4,6,c(1,2,3),-0.2,-0.3), 3)
  expect_length(object = qkappa(c(.2,.3,.6),6,c(1,2,3,1,2,3),-1,-0.3), 6)
  expect_length(object = rkappa(n = 15,6,c(0.01,3),-0.2,-0.2), 15)
  ## test added for GEV and GLO for sensible values in pdf and cdf when x-values are out of the domain 
  expect_lte(pgev(2*(40+5/0.2),40,5,0.2),1) ## cdf maximum value is 1 for points after upper bound
  expect_gte(pgev((40-5/0.2)-20,40,5,-0.2),0) ## cdf has value 0 for points below lower bound 
  expect_lte(dgev(2*(40+5/0.2),40,5,0.2),10^-7) ## cdf maximum value is 1 for points after upper bound
  expect_lte(dgev((40-5/0.2)-10,40,5,-0.2),10^-7) ## cdf has value 0 for points below lower bound 
  expect_lte(pglo(2*(40+5/0.2),40,5,0.2),1) ## cdf maximum value is 1 for points after upper bound
  expect_gte(pglo((40-5/0.2)-20,40,5,-0.2),0) ## cdf has value 0 for points below lower bound 
  expect_lte(dglo(2*(40+5/0.2),40,5,0.2),10^-7) ## cdf maximum value is 1 for points after upper bound
  expect_lte(dglo((40-5/0.2)-10,40,5,-0.2),10^-7) ## cdf has value 0 for points below lower bound 
  ## check that quantiles are as those in lmom
  expect_equal(qkappa(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2, -0.4), 
               lmom::quakap(c(0.5,0.99,0.995,0.995,0.999),c(4,6,0.2, -0.4)))
  expect_equal(qkappa(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2, -0.4, lower.tail = FALSE), 
               lmom::quakap(1-c(0.5,0.99,0.995,0.995,0.999),c(4,6,0.2, -0.4)))
  expect_equal(qglo(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2), 
               lmom::quaglo(c(0.5,0.99,0.995,0.995,0.999),c(4,6,0.2)))
  expect_equal(qglo(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2, lower.tail = FALSE), 
               lmom::quaglo(1-c(0.5,0.99,0.995,0.995,0.999),c(4,6,0.2)))
  # notice the GEV has the same parametrisation as in Hosking and Wallis
  expect_equal(qgev(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2), 
               lmom::quagev(c(0.5,0.99,0.995,0.995,0.999),c(4,6,0.2)))
  expect_equal(qgev(c(0.5,0.99,0.995,0.995,0.999),4,6,0.2, lower.tail = FALSE), 
               lmom::quagev(1-c(0.5,0.99,0.995,0.995,0.999),c(4,6,0.2)))
})

set.seed(124); y <- rglo(n = 600, 60, 5, 0.2)
xx <- matrix(runif(600), ncol=1)

test_that("gev fitting works", {
  fit1 <- gevd.fit(y, show = FALSE)
  expect_length(object = fit1$mle, 3)
  fit2 <- gevd.fit(y, show = FALSE, fixedPars = list(sh=-0.4))
  expect_length(object = fit2$mle, 2)
  fit3 <- gevd.fit(y, show = FALSE, 
                   ydat = xx, mul=1, 
                   fixedPars = list(sh=-0.4))
  # fit3$mle
  # [1] 58.007815 -1.183047  5.776229
  expect_length(object = fit3$mle, 3)
  expect_gt(fit3$mle[1], 56.671)
  expect_lt(fit3$mle[1], 56.672)
})



test_that("glo fitting works", {
  fit1 <- glod.fit(y, show = FALSE)
  expect_length(object = fit1$mle, 3)
  expect_gt(fit1$mle[3], 0.183)
  expect_lt(fit1$mle[3], 0.184)
  fit2 <- glod.fit(y, show = FALSE, fixedPars = list(sh=0.2))
  expect_length(object = fit2$mle, 2)
  fit3 <- glod.fit(y, show = FALSE, 
                   ydat = xx, mul=1, 
                   fixedPars = list(sh=0.2))
  # fit3$mle
  # [1] 61.170951 -1.720171  4.765090
  expect_length(object = fit3$mle, 3)
  expect_gt(fit3$mle[1], 61.170)
  expect_lt(fit3$mle[1], 61.171)
})



test_that("glo - cv model- fitting works", {
  fit1 <- glocvd.fit(y, show = FALSE)
  expect_length(object = fit1$mle, 3)
  expect_gt(fit1$mle[3], 0.183)
  expect_lt(fit1$mle[3], 0.184)
  fit2 <- glocvd.fit(y, show = FALSE, fixedPars = list(sh=0.2))
  expect_length(object = fit2$mle, 2)
  fit3 <- glocvd.fit(y, show = FALSE, 
                   ydat = xx, mul=1, 
                   fixedPars = list(sh=0.2))
  # fit3$mle
  # [1] 61.1824128 -1.7442033  0.0789351
  expect_length(object = fit3$mle, 3)
  expect_gt(fit3$mle[1], 61.1824)
  expect_lt(fit3$mle[1], 61.1825)
})

test_that("kappa fitting works", {
  fit1 <- kappad.fit(y, show = FALSE)
  expect_length(object = fit1$mle, 4)
  expect_gt(fit1$mle[3], 0.154)
  expect_lt(fit1$mle[3], 0.155)
  fit2 <- kappad.fit(y, show = FALSE, fixedPars = list(sh = 0.2, sh2=-1))
  expect_length(object = fit2$mle, 2)
  fit3 <- kappad.fit(y, show = FALSE, 
                   ydat = xx, mul=1, 
                   fixedPars = list(sh=0.2, sh2 = -1))
  # fit3$mle
  # [1]  61.170951 -1.720171  4.765090
  expect_length(object = fit3$mle, 3)
  expect_gt(fit3$mle[1], 61.170)
  expect_lt(fit3$mle[1], 61.171)
})

