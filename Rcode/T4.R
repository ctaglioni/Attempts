library(methods)
#library(nzreg)
library(dplyr)
library(docopt)
library(demest)
library(demItaly)
library(abind)
library(beepr)

dimnames(italy.popn.reg)$time <- 2005:2016
census_year_erp<- Counts(italy.popn.reg, dimscales = c(time="Points")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = "time")

#dim(italy.popn.reg)
#dim(census_year_erp)

population <- census_year_erp

reg_births<- Counts(italy.births.reg, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = "time")

reg_deaths <- Counts(italy.deaths.reg, dimscales = c(time="Intervals")) %>%
  collapseDimension(margin = "time")

arrivals <- Counts(italy.ext.in , dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = "time")

departures <- Counts(italy.ext.out, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = "time")



account<- Movements(population = population,
                    births = reg_births,
                    entries = list(external_in = arrivals),
                    exits = list(deaths = reg_deaths,
                                 external_out = departures)) %>%
  makeConsistent()




systemModels <- list(Model(population ~ Poisson(mean ~ time, useExpose = FALSE),
                           time ~ DLM(level = Level(scale = HalfT(scale = 0.05)),
                                      trend = Trend(scale = HalfT(scale = 0.05)),
                                      damp = NULL,
                                      error = Error(scale = HalfT(scale = 0.05))),
                           jump = 0.0005),
                     Model(births ~ Poisson(mean ~ time),
                           time ~ DLM(level = Level(scale = HalfT(scale = 0.05)),
                                      trend = NULL,
                                      damp = NULL,
                                      error = Error(scale = HalfT(scale = 0.05))),
                           jump = 0.005),
                     Model(deaths ~ Poisson(mean ~ time),
                           time ~ DLM(trend = NULL, damp = NULL),
                           jump = 0.005),
                     Model(external_in ~ Poisson(mean ~ time),
                           time ~ DLM(level = Level(scale = HalfT(scale = 0.3)),
                                      trend = NULL,
                                      damp = NULL,
                                      error = Error(scale = HalfT(scale = 0.3))),
                           jump = 0.005),
                     Model(external_out ~ Poisson(mean ~ time),
                           time ~ DLM(level = Level(scale = HalfT(scale = 0.3)),
                                      trend = Trend(scale = HalfT(scale = 0.3)),
                                      damp = NULL,
                                      error = Error(scale = HalfT(scale = 0.3))),
                           jump = 0.005))

datasets <- list(census_year_erp = census_year_erp,
                 reg_births = reg_births,
                 reg_deaths = reg_deaths,
                 arrivals = arrivals,
                 departures = departures)

sd <- census_year_erp %>%
  collapseDimension(margin = "time") %>%
  as("Values") 
sd<- sd * 0.0025
mean <- sd / sd


dataModels <- list(Model(census_year_erp ~ NormalFixed(mean = mean, sd = sd),
                         series = "population"),
                   Model(reg_births ~ PoissonBinomial(prob = 0.98),
                         series = "births"),
                   Model(reg_deaths ~ PoissonBinomial(prob = 0.95),
                         series = "deaths"),
                   Model(arrivals ~ Poisson(mean ~ 1),
                         `(Intercept)` ~ ExchFixed(sd = 0.05),
                         lower = 0.8,
                         upper = 1.2,
                         priorSD = HalfT(scale = 0.075),
                         series = "external_in",
                         jump = 0.005),
                   Model(departures ~ Poisson(mean ~ 1),
                         `(Intercept)` ~ ExchFixed(sd = 0.025),
                         lower = 0.8,
                         upper = 1.2,
                         priorSD = HalfT(scale = 0.05),
                         series = "external_out",
                         jump = 0.005))

filename <- "C:/0_PhD/Thesis/Thesis_R/onedim1.est"

n_sim <- 200000
n_burnin <- 200000
n_chain <- 4
n_thin <- 500



beep(sound = 2, estimateAccount(account = account,
                                systemModels = systemModels,
                                datasets = datasets,
                                dataModels = dataModels,
                                filename = filename,
                                nBurnin = n_burnin,
                                nSim = n_sim,
                                nChain = n_chain,
                                nThin = n_thin,
                                useC = TRUE)
)



fetchSummary(filename)

pop.chain<- fetch(filename, where=c("account", "population"))
bir.chain<- fetch(filename, where=c("account", "births"))
dea.chain <- fetch(filename, where=c("account", "deaths"))
imm.chain <- fetch(filename, where=c("account", "external_in"))
emi.chain <- fetch(filename, where=c("account", "external_out"))

firstrow<- 0.5*(n_sim/n_thin*n_chain)
lastrow<- dim(pop.chain)[2]

dim(pop.chain)
popest<- apply(pop.chain[,firstrow:lastrow], 1, mean)
birest <- apply(bir.chain[,firstrow:lastrow], 1, mean)
deaest <- apply(dea.chain[,firstrow:lastrow], 1, mean)
immest <- apply(imm.chain[,firstrow:lastrow], 1, mean)
emiest <- apply(emi.chain[,firstrow:lastrow], 1, mean)


#-----------------------------------
dplot( ~ time, data = popest - population,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Population differences 2006-2015")

plot(popest - population, ylab="Estimate - Data", 
     main = "Population: estimate vs data")

dplot( ~ time, data = birest - reg_births,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births differences 2006-2015")
plot(birest- reg_births, ylab="Estimate - Data", 
     main = "Births: estimate vs data")

dplot( ~ time, data = deaest - reg_deaths,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Deaths differences 2006-2015")
plot(deaest - reg_deaths, ylab="Estimate - Data", 
     main = "Deaths: estimate vs data")

dplot( ~ time, data = immest - arrivals,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Immigration differences 2006-2015")
plot(immest - arrivals, ylab="Estimate - Data", 
     main = "Immigration: estimate vs data")

dplot( ~ time, data = emiest - departures,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Emigration differences 2006-2015")
plot(emiest - departures, ylab="Estimate - Data", 
     main = "Emigration: estimate vs data")


#-----------------------------------

dplot( ~ time, data = pop.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Population estimation 2006-2015",
       overlay = list(values = population, col = "red", lwd= 2))

dplot( ~ time, data = bir.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation 2006-2015",
       overlay = list(values = reg_births, col = "red", lwd= 2))

dplot( ~ time, data = dea.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Deaths estimation 2006-2015",
       overlay = list(values = reg_deaths, col = "red", lwd= 2))

dplot( ~ time, data = imm.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Immigration estimation 2006-2015",
       overlay = list(values = arrivals, col = "red", lwd= 2))

dplot( ~ time, data = emi.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Emigration estimation 2006-2015",
       overlay = list(values = departures, col = "red", lwd= 2))

population
#------------------------------

MCMC <- fetchMCMC(filename)
names(MCMC)
library(coda)

plot(MCMC$account.population, sub=names(MCMC)[1])
plot(MCMC$account.births, sub=names(MCMC)[2])
plot(MCMC$account.external_in, sub=names(MCMC)[3]) # problems in 2014
plot(MCMC$account.deaths, sub=names(MCMC)[4])
plot(MCMC$account.external_out, sub=names(MCMC)[5]) # problems especially in 2014
plot(MCMC$systemModels.population.likelihood.count, sub=names(MCMC)[6])
plot(MCMC$systemModels.population.prior.mean, sub=names(MCMC)[7])
plot(MCMC$systemModels.population.hyper.time.scaleLevel, sub=names(MCMC)[8])
plot(MCMC$systemModels.population.hyper.time.scaleTrend, sub=names(MCMC)[9])
plot(MCMC$systemModels.population.hyper.time.scaleError, sub=names(MCMC)[10])
plot(MCMC$systemModels.births.likelihood.rate, sub=names(MCMC)[11])
plot(MCMC$systemModels.births.prior.mean, sub=names(MCMC)[12])
plot(MCMC$systemModels.births.hyper.time.scaleLevel, sub=names(MCMC)[13])
plot(MCMC$systemModels.births.hyper.time.scaleError, sub=names(MCMC)[14])
plot(MCMC$systemModels.external_in.likelihood.rate, sub=names(MCMC)[15]) #2014! But seems to converge at the end
plot(MCMC$systemModels.external_in.prior.mean, sub=names(MCMC)[16]) # prior okay! Converging.
plot(MCMC$systemModels.external_in.hyper.time.scaleLevel, sub=names(MCMC)[17])
plot(MCMC$systemModels.external_in.hyper.time.scaleError, sub=names(MCMC)[18])
plot(MCMC$systemModels.deaths.likelihood.rate, sub=names(MCMC)[19])
plot(MCMC$systemModels.deaths.prior.mean, sub=names(MCMC)[20])
plot(MCMC$systemModels.deaths.hyper.time.scaleLevel, sub=names(MCMC)[21])
#plot(MCMC$systemModels.deaths.hyper.time.scaleTrend, sub=names(MCMC)[22])
plot(MCMC$systemModels.deaths.hyper.time.scaleError, sub=names(MCMC)[22])
plot(MCMC$systemModels.external_out.likelihood.rate, sub=names(MCMC)[23])
plot(MCMC$systemModels.external_out.prior.mean, sub=names(MCMC)[24]) # 2014! But seems to converge at the end
plot(MCMC$systemModels.external_out.hyper.time.scaleLevel, sub=names(MCMC)[25]) # prior ok
# plot(MCMC$systemModels.external_out.hyper.time.scaleTrend, sub=names(MCMC)[27])
plot(MCMC$systemModels.external_out.hyper.time.scaleError, sub=names(MCMC)[26])
plot(MCMC$dataModels.arrivals.likelihood.rate, sub=names(MCMC)[27]) # as system.model
plot(MCMC$dataModels.arrivals.prior.mean, sub=names(MCMC)[28])
plot(MCMC$dataModels.arrivals.prior.sd, sub=names(MCMC)[29]) #...convergence?
plot(MCMC$dataModels.departures.likelihood.rate, sub=names(MCMC)[30]) # most problematic chain
plot(MCMC$dataModels.departures.prior.mean, sub=names(MCMC)[31]) # prior ok
plot(MCMC$dataModels.departures.prior.sd, sub=names(MCMC)[32])

# plot(MCMC$dataModels.departures.hyper.time.scaleLevel,sub=names(MCMC)[34])
# plot(MCMC$dataModels.departures.hyper.time.scaleTrend, sub=names(MCMC)[35])
# plot(MCMC$dataModels.departures.hyper.time.scaleError, sub=names(MCMC)[36])

