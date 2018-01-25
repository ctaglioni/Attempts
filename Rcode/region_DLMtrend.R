#-------------------
# Only region
#-------------------

library(methods)
#library(nzreg)
library(dplyr)
library(docopt)
library(demest)
library(demItaly)

setwd("C:\\Users\\Cha\\Documents\\demItaly\\data")

census_year_erp<- Counts(italy.popn.reg, dimscales = c(time="Points")) %>%
  collapseDimension(margin = c("time", "region"))%>%
  subarray(time != c("2016", "2017"))

population<- extrapolate(census_year_erp, labels = "2005")


reg_births<- Counts(italy.births.reg, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = c("time", "region"))
plot(reg_births)

reg_deaths <- Counts(italy.deaths.reg, dimscales = c(time="Intervals")) %>%
  collapseDimension(margin = c("time", "region"))
plot(reg_deaths)

arrivals<- Counts(italy.ext.in , dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = c("time", "region"))

departures<- Counts(italy.ext.out, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = c("time", "region"))

intarr<- Counts(italy.int.in , dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = c("time", "region"))

intdep<- Counts(italy.int.out, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = c("time", "region"))



# Regional means
Regions <- dimnames(population)$region
mean_pop <- data.frame(region = Regions, mean = apply(population, 2, mean))
mean_bir <- data.frame(region = Regions, mean=apply(reg_births, 2, mean))
mean_dea <- data.frame(region = Regions, mean=apply(reg_deaths, 2, mean))
mean_arr <- data.frame(region = Regions, mean=apply(arrivals, 2, mean))
mean_dep <- data.frame(region = Regions, mean=apply(departures, 2, mean))
mean_intarr <- data.frame(region = Regions, mean=apply(intarr, 2, mean))
mean_intdep <- data.frame(region = Regions, mean=apply(intdep, 2, mean))

account<- Movements(population = population,
                    births = reg_births,
                    entries = list(external_in = arrivals,
                                   internal_in = intarr),
                    exits = list(deaths = reg_deaths,
                                 external_out = departures,
                                 internal_out = intdep)) %>%
  makeConsistent()

systemModels<- list(Model(population ~ Poisson(mean ~ time + region, useExpose=FALSE),
                          time ~ DLM(trend = Trend(initial = Initial(sd = 2),
                                           scale = HalfT(df = 10, max = 2))),
                          region ~ Exch(),
                          jump = 0.003),
                    Model(births ~ Poisson(mean ~ time + region),
                          time ~ DLM(trend = Trend(initial = Initial(sd = 2),
                                                   scale = HalfT())),
                          region ~ Exch(),
                          jump = 0.02),
                    Model(deaths ~ Poisson(mean ~ time + region),
                          time ~ DLM(trend = Trend(initial = Initial(sd = 2),
                                                   scale = HalfT())),
                          region ~ Exch(),
                          jump = 0.015),
                    Model(external_in ~ Poisson(mean ~ time + region),
                          time ~ DLM(trend = Trend(initial = Initial(sd = 2),
                                                           scale = HalfT())), 
                          region ~ Exch(error = Error(robust = TRUE)),
                          jump = 0.03),
                    Model(external_out ~ Poisson(mean ~ time + region),
                          time ~ DLM(trend = Trend(initial = Initial(sd = 2),
                                                           scale = HalfT())),
                          region ~ Exch(error = Error(robust = TRUE)),
                          jump = 0.03),
                    Model(internal_in ~ Poisson(mean ~ time + region),
                          time ~ DLM(trend = Trend(initial = Initial(sd = 2),
                                                           scale = HalfT())), 
                          region ~ Exch(error = Error(robust = TRUE)),
                          jump = 0.03),
                    Model(internal_out ~ Poisson(mean ~ time + region),
                          time ~ DLM(trend = Trend(initial = Initial(sd = 2),
                                                           scale = HalfT())),
                          region ~ Exch(error = Error(robust = TRUE)),
                          jump = 0.03))

# the population size of the region is not directly linked to the number of migrants
# i.e. only exchangeable prior, nor covariate

datasets <- list(census_year_erp = census_year_erp,
                 reg_births = reg_births,
                 reg_deaths = reg_deaths,
                 arrivals = arrivals,
                 departures = departures,
                 intarr = intarr,
                 intdep = intdep)

dataModels <- list(Model(census_year_erp ~ PoissonBinomial(prob = 0.95),
                         series = "population"),
                   Model(reg_births ~ PoissonBinomial(prob = 0.98),
                         series = "births"),
                   Model(reg_deaths ~ PoissonBinomial(prob = 0.98),
                         series = "deaths"),
                   Model(arrivals ~ Poisson(mean ~ 1),
                         series = "external_in",
                         jump = 0.05),
                   Model(departures ~ Poisson(mean ~ 1),
                         series = "external_out",
                         jump = 0.05),
                   Model(intarr ~ Poisson(mean ~ 1),
                         series = "internal_in",
                         jump = 0.05),
                   Model(intdep ~ Poisson(mean ~ 1),
                         series = "internal_out",
                         jump = 0.05))

onlyreg <- "C:/0_PhD/Thesis/Thesis_R/RegionDLM1.est"

n_sim <- 10000
n_burnin <- 10000
n_chain <- 4
n_thin <- 50

estimateAccount(account = account,
                systemModels = systemModels,
                datasets = datasets,
                dataModels = dataModels,
                filename = onlyreg,
                nBurnin = n_burnin,
                nSim = n_sim,
                nChain = n_chain,
                nThin = n_thin,
                useC = TRUE)

fetchSummary(onlyreg)
showModel(onlyreg)
listContents(onlyreg)

pop.chain<- fetch(onlyreg, where=c("account", "population"))
bir.chain <- fetch(onlyreg, where=c("account", "births"))
dea.chain <- fetch(onlyreg, where=c("account", "deaths"))
imm.chain <- fetch(onlyreg, where=c("account", "external_in"))
emi.chain <- fetch(onlyreg, where=c("account", "external_out"))
iin.chain <- fetch(onlyreg, where=c("account", "internal_in"))
ein.chain <- fetch(onlyreg, where=c("account", "internal_out"))

firstrow <- 0.5*(n_sim/n_thin*n_chain)
lastrow <- dim(pop.chain)[3]

dim(pop.chain[,,firstrow:lastrow])
popest <- apply(pop.chain[,,firstrow:lastrow], c(1,2), mean)
birest <- apply(bir.chain[,,firstrow:lastrow], c(1,2), mean)
deaest <- apply(dea.chain[,,firstrow:lastrow], c(1,2), mean)
immest <- apply(imm.chain[,,firstrow:lastrow], c(1,2), mean)
emiest <- apply(emi.chain[,,firstrow:lastrow], c(1,2), mean)
iinest <- apply(iin.chain[,,firstrow:lastrow], c(1,2), mean)
einest <- apply(ein.chain[,,firstrow:lastrow], c(1,2), mean)

# Estimate - truth
dplot( ~ time | region, data = popest - population,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Population differences 2006-2015")
dplot( ~ time | region, data = birest - reg_births,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births differences 2006-2015")
dplot( ~ time | region, data = deaest - reg_deaths,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Deaths differences 2006-2015")
dplot( ~ time | region, data = immest - arrivals,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Immigration differences 2006-2015")
dplot( ~ time | region, data = emiest - departures,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Emigration differences 2006-2015")
dplot( ~ time | region, data = iinest - intarr,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal imm differences 2006-2015")
dplot( ~ time | region, data = einest - intdep,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal emi differences 2006-2015")

# Estimation by year
dplot( ~ time | region, data = pop.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Population estimation 2006-2014",
       overlay = list(values = population, col = "red", lwd= 2))

dplot( ~ time | region, data = bir.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation 2006-2014",
       overlay = list(values = reg_births, col = "red", lwd= 2))

dplot( ~ time | region, data = dea.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Deaths estimation 2006-2014",
       overlay = list(values = reg_deaths, col = "red", lwd= 2))

dplot( ~ time | region, data = imm.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Immigration estimation 2006-2014",
       overlay = list(values = arrivals, col = "red", lwd= 2))

dplot( ~ time | region, data = emi.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Emigration estimation 2006-2014",
       overlay = list(values = departures, col = "red", lwd= 2))

dplot( ~ time | region, data = iin.chain ,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal immigration estimation 2006-2014",
       overlay = list(values = intarr, col = "red", lwd= 2))

dplot( ~ time | region, data = ein.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal emigration estimation 2006-2014",
       overlay = list(values = intdep, col = "red", lwd= 2))


MCMC <- fetchMCMC(onlyreg)
names(MCMC)
library(coda)
plot(MCMC$account.population)
plot(MCMC$systemModels.births.likelihood.rate)
gelman.diag(MCMC$systemModels.population.likelihood.count)
