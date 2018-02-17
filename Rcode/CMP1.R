library(methods)
#library(nzreg)
library(dplyr)
library(docopt)
library(demest)
library(demItaly)
library(abind)
library(beepr)

y <- Counts(italy.deaths.reg, dimscales = c(time="Intervals")) %>%
  collapseDimension(margin = "time")

#dimnames(italy.popn.reg)$time <- 2001:2016
expose<- Counts(italy.popn.reg, dimscales = c(time="Intervals")) %>%
  subarray(time != "2016") %>%
  subarray(time > "2005") %>%
  collapseDimension(margin = "time")

filename <- tempfile()
estimateModel(Model(y ~ CMP(mean ~ 1,
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 1),
                                                    scale = HalfT(scale = 1)),
                            useExpose = FALSE),
                    jump = 0.001),
              y = y,
              
              filename = filename,
              nBurnin = 0,
              nSim = 200,
              nChain = 4,
              nThin = 1)

fetchSummary(filename)


MCMC <-fetchMCMC(filename)
plot(MCMC$model.likelihood.count)
plot(MCMC$model.likelihood.dispersion )
plot(MCMC$model.prior.dispersion.mean)
plot(MCMC$model.prior.dispersion.sd )
plot(MCMC$model.hyper.age.scaleLevel )
plot(MCMC$model.prior.count.mean )
plot(MCMC$model.prior.count.sd )
plot(MCMC$model.hyper.age.scaleTrend)




#---------------------------------
# Model after correction
#---------------------------------

library(dplyr)
library(demest)

y <- demdata::uk.deaths %>%
  Counts(dimscales = c(year = "Intervals")) %>%
  subarray(year > 1950) %>%
  collapseIntervals(dimension = "age", breaks = seq(0, 90, 5))

expose <- demdata::uk.exposure %>%
  Counts(dimscales = c(year = "Intervals")) %>%
  subarray(year > 1950) %>%
  collapseIntervals(dimension = "age", breaks = seq(0, 90, 5))
filename <- tempfile()
estimateModel(Model(y ~ CMP(mean ~ age + sex + year,
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 0.0001),
                                                    scale = HalfT(scale = 0.0001))),
                    age ~ DLM(damp = NULL),
                    jump = 0.000005),
              y = y,
              exposure = expose,
              filename = filename,
              nBurnin = 2000,
              nSim = 2000,
              nChain = 4,
              nThin = 20)
fetchSummary(filename)

MCMC <-fetchMCMC(filename)
plot(MCMC$model.likelihood.rate)
plot(MCMC$model.likelihood.dispersion )
plot(MCMC$model.prior.dispersion.mean)
plot(MCMC$model.prior.dispersion.sd )
plot(MCMC$model.hyper.age.scaleLevel )
plot(MCMC$model.prior.rate.mean )
plot(MCMC$model.prior.rate.sd )
plot(MCMC$model.hyper.age.scaleTrend)
plot(MCMC$model.hyper.age.scaleError )
plot(MCMC$model.hyper.year.damp )
plot(MCMC$model.hyper.year.scaleError )


#----------------------------
# Old checking
#----------------------------
library(dplyr)
library(demest)
library(beepr)

y <- demdata::uk.deaths %>%
  Counts(dimscales = c(year = "Intervals")) %>%
  subarray(year > 1950) %>%
  collapseIntervals(dimension = "age", breaks = seq(0, 90, 5))
expose <- demdata::uk.exposure %>%
  Counts(dimscales = c(year = "Intervals")) %>%
  subarray(year > 1950) %>%
  collapseIntervals(dimension = "age", breaks = seq(0, 90, 5))

filename <- tempfile()
estimateModel(Model(y ~ CMP(mean ~ age + sex + year,
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 0.1),
                                                    scale = HalfT(scale = 0.001))),
                    age ~ DLM(damp = NULL),
                    upper = 2,
                    jump = 0.1),
              y = y,
              exposure = expose,
              filename = filename,
              nBurnin = 20000,
              nSim = 20000,
              nChain = 4,
              nThin = 200)
fetchSummary(filename)
showModel(filename)
MCMC <-fetchMCMC(filename)
plot(MCMC$model.likelihood.rate)
plot(MCMC$model.likelihood.dispersion )
plot(MCMC$model.prior.dispersion.mean)
plot(MCMC$model.prior.dispersion.sd )
plot(MCMC$model.hyper.age.scaleLevel )
plot(MCMC$model.prior.rate.mean )
plot(MCMC$model.prior.rate.sd )
plot(MCMC$model.hyper.age.scaleTrend)
plot(MCMC$model.hyper.age.scaleError )
plot(MCMC$model.hyper.year.damp )
plot(MCMC$model.hyper.year.scaleError )


#------------------------------------------
y <- demdata::uk.deaths %>%
  Counts(dimscales = c(year = "Intervals")) %>%
  subarray(year > 1950) %>%
  collapseDimension(margin = "year")

expose <- demdata::uk.exposure %>%
  Counts(dimscales = c(year = "Intervals")) %>%
  subarray(year > 1950) %>%
  collapseDimension(margin = "year")


beep(sound = 2, estimateModel(Model(y ~ CMP(mean ~ year,
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 1),
                                                    scale = HalfT(scale = 0.1))),
                    year ~ DLM(),
                    upper = 2,
                    jump = 0.0000001),
              y = y,
              exposure = expose,
              filename = filename,
              nBurnin = 20000,
              nSim = 20000,
              nChain = 4,
              nThin = 200))

beep(sound = 2, estimateModel(Model(y ~ CMP(mean ~ year,
                                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 0.01),
                                                                    scale = HalfT(scale = 0.001))),
                                    year ~ DLM(),
                                    upper = 2,
                                    jump = 0.0000001),
                              y = y,
                              exposure = expose,
                              filename = filename,
                              nBurnin = 2000000,
                              nSim = 2000000,
                              nChain = 4,
                              nThin = 2000))


fetchSummary(filename)
showModel(filename)
listContents(filename)
fetch(filename, c("model","likelihood","noProposalY"))
MCMC <-fetchMCMC(filename)
plot(MCMC$model.likelihood.rate)
plot(MCMC$model.likelihood.dispersion )
plot(MCMC$model.prior.dispersion.mean)
plot(MCMC$model.prior.dispersion.sd )
plot(MCMC$model.prior.rate.mean )
plot(MCMC$model.hyper.year.damp )
plot(MCMC$model.hyper.year.scaleError )
plot(MCMC$model.hyper.year.scaleLevel)
plot(MCMC$model.hyper.year.scaleTrend)
