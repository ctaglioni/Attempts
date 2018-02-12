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
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 0.1),
                                                    scale = HalfT(scale = 0.001))),
                    age ~ DLM(damp = NULL),
                    upper = 2,
                    jump = 0.1),
              y = y,
              exposure = expose,
              filename = filename,
              nBurnin = 200000,
              nSim = 200000,
              nChain = 4,
              nThin = 500)
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

