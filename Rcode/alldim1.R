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
MigrAgeGroupDisag <- function(migr, s = 2, r = 22, y = 11){ # migration dataset
  # age group width according ISTAT classification
  ag.width <- c(18, 22, 25, 26) # "65+" has width 26 i.e. "65-90"
  # average per single year age group
  single <- apply(migr, c(2,3,4), "/", ag.width)
  # 5y age groups reconstitution
  ag014<- floor(single[1,,,]*5)
  ag1519 <- floor(single[1,,,]*3+single[2,,,]*2)
  ag2039 <- floor(single[2,,,]*5)
  ag4064 <- floor(single[3,,,]*5)
  ag6589 <- floor(single[4,,,]*5)
  # last age group "90+" includes what is left from the others
  ag90 <- colSums(migr) - 3*ag014 - ag1519 - 4*ag2039 - 5*ag4064 - 5*ag6589
  # new array
  arr <- array(NA, dim = c(19,s,r,y))
  
  a1<- aperm(abind(ag014,ag014,ag014, along = 4), c(4,1,2,3))
  a2<- aperm(abind(ag1519, along = 4), c(4,1,2,3))
  a3<- aperm(abind(ag2039,ag2039,ag2039,ag2039, along = 4), c(4,1,2,3))
  a4<- aperm(abind(ag4064,ag4064,ag4064,ag4064,ag4064, along = 4), c(4,1,2,3))
  a5<- aperm(abind(ag6589,ag6589,ag6589, ag6589, ag6589, along = 4), c(4,1,2,3))
  a6<- aperm(abind(ag90,along = 4), c(4,1,2,3))
  arr <- abind(a1,a2,a3,a4,a5,a6, along = 1)
  return(arr)
}

census_year_erp<- Counts(italy.popn.reg, dimscales = c(time="Points")) %>%
  subarray(time != c("2016", "2017"))%>%
  collapseDimension(margin = c("age", "region", "time"))

population2<- extrapolate(census_year_erp, labels = "2005")
population<- population2[ , , c(1,6,11)]
dim(population)

reg_births<- Counts(italy.births.reg, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016")%>%
  reallocateToEndAges() %>%
  collapseIntervals( dimension= "time", width = 5) %>%
  collapseIntervals(dimension = "age", width = 5)
plot(reg_births)

reg_deaths<- Counts(italy.deaths.reg, dimscales = c(time="Intervals"))%>%
  collapseDimension(margin = c("age", "region", "time")) %>%
  collapseIntervals(dimension = "time", width = 5)
plot(reg_deaths)

arr <- MigrAgeGroupDisag(italy.ext.in)
dim(arr)
dimnames(arr) <-list(age=dimnames(italy.popn.reg)[[1]],
                     sex=dimnames(italy.popn.reg)[[2]],
                     region = dimnames(italy.popn.reg)[[3]],
                     time = c(2006:2016))

arrivals<- Counts(arr , dimscales = c(time="Intervals")) %>%
  subarray(time !="2016")%>%
  collapseDimension(margin = c("age", "region", "time")) %>%
  collapseIntervals(dimension = "time", width = 5)
plot(arrivals)

dep <-MigrAgeGroupDisag(italy.ext.out)
dim(dep)
dimnames(dep) <-list(age=dimnames(italy.popn.reg)[[1]],
                     sex=dimnames(italy.popn.reg)[[2]],
                     region = dimnames(italy.popn.reg)[[3]],
                     time = c(2006:2016))

departures<- Counts(dep, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016")%>%
  collapseDimension(margin = c("age", "region", "time")) %>%
  collapseIntervals(dimension = "time", width = 5)
plot(departures)

iin <-MigrAgeGroupDisag(italy.int.in)
dim(iin)
dimnames(iin) <-list(age=dimnames(italy.popn.reg)[[1]],
                     sex=dimnames(italy.popn.reg)[[2]],
                     region = dimnames(italy.popn.reg)[[3]],
                     time = c(2006:2016))

intarr<- Counts(iin , dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = c("age", "region", "time")) %>%
  collapseIntervals(dimension = "time", width = 5)

ein <-MigrAgeGroupDisag(italy.int.out)
dim(ein)
dimnames(ein) <-list(age=dimnames(italy.popn.reg)[[1]],
                     sex=dimnames(italy.popn.reg)[[2]],
                     region = dimnames(italy.popn.reg)[[3]],
                     time = c(2006:2016))


intdep<- Counts(ein, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016")%>%
  collapseDimension(margin = c("age", "region", "time")) %>%
  collapseIntervals(dimension = "time", width = 5)



account<- Movements(population = population,
                    births = reg_births,
                    entries = list(external_in = arrivals,
                                   internal_in = intarr),
                    exits = list(deaths = reg_deaths,
                                 external_out = departures,
                                 internal_out = intdep)) %>%
  makeConsistent()

systemModels<- list(Model(population ~ Poisson(mean ~ age + time + region, useExpose=FALSE),
                          time ~ DLM(damp = NULL),
                          region ~ Exch(),
                          age ~ DLM(),
                          jump = 0.003),
                    Model(births ~ Poisson(mean ~ age + time + region),
                          time ~ DLM(trend = NULL, damp = NULL),
                          region ~ Exch(),
                          age ~ DLM(),
                          jump = 0.02),
                    Model(deaths ~ Poisson(mean ~ age + time + region),
                          time ~ DLM(trend = NULL, damp = NULL),
                          region ~ Exch(),
                          age ~ DLM(),
                          jump = 0.015),
                    Model(external_in ~ Poisson(mean ~ age + time + region),
                          time ~ DLM(trend = NULL), 
                          region ~ Exch(),
                          age ~ DLM(),
                          jump = 0.03),
                    Model(external_out ~ Poisson(mean ~ age + time + region),
                          time ~ DLM(trend = NULL),
                          region ~ Exch(),
                          age ~ DLM(),
                          jump = 0.03),
                    Model(internal_in ~ Poisson(mean ~ age + time + region),
                          time ~ DLM(trend = NULL), 
                          region ~ Exch(),
                          age ~ DLM(),
                          jump = 0.03),
                    Model(internal_out ~ Poisson(mean ~ age + time + region),
                          time ~ DLM(trend = NULL),
                          region ~ Exch(),
                          age ~ DLM(),
                          jump = 0.03))

# the population size of the region is not directly linked to the number of migrants
# i.e. only exchangeable prior, nor covariate

datasets <- list(census_year_erp = population,
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

onlyreg <- "C:/0_PhD/Thesis/Thesis_R/alldim1.est"

n_sim <- 1000
n_burnin <- 1000
n_chain <- 4
n_thin <- 10

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
lastrow <- dim(pop.chain)[length(dim(pop.chain))]

dim(pop.chain[,,,firstrow:lastrow])
popest <- apply(pop.chain[,,,firstrow:lastrow], c(1,2,3), mean)
birest2<- apply(bir.chain[,,,,firstrow:lastrow], c(1,2,3,4), mean)
birest<- apply(birest2, c(1,2,3),sum)
deaest2 <- apply(dea.chain[,,,,firstrow:lastrow], c(1,2,3,4), mean)
deaest<- apply(deaest2, c(1,2,3),sum)
immest2 <- apply(imm.chain[,,,,firstrow:lastrow], c(1,2,3,4), mean)
immest<- apply(immest2, c(1,2,3),sum)
emiest2 <- apply(emi.chain[,,,,firstrow:lastrow], c(1,2,3,4), mean)
emiest<- apply(emiest2, c(1,2,3),sum)
iinest2 <- apply(iin.chain[,,,,firstrow:lastrow], c(1,2,3,4), mean)
iinest<- apply(iinest2, c(1,2,3),sum)
einest2 <- apply(ein.chain[,,,,firstrow:lastrow], c(1,2,3,4), mean)
einest<- apply(einest2, c(1,2,3),sum)


# Estimate - truth
dplot( ~ age | region, data = popest - population,
       subarray=time=="2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Population differences 2015")
dplot( ~ age | region, data = birest - aperm(reg_births, c(2,1,3)),
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births differences 2015")
dplot( ~ age | region, data = deaest - reg_deaths,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Deaths differences 2015")
dplot( ~ age | region, data = immest - arrivals,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Immigration differences 2015")
dplot( ~ age | region, data = emiest - departures,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Emigration differences 2015")
dplot( ~ age | region, data = iinest - intarr,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal imm differences 2015")
dplot( ~ age | region, data = einest - intdep,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal emi differences 2015")


# Estimation by year
dplot( ~ age | region, data = pop.chain,
       subarray=time=="2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Population estimation 2006-2014",
       overlay = list(values = subarray(population, time == "2015"), col = "red", lwd= 2))

dplot( ~ age | region, data = bir.chain,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation 2006-2014",
       overlay = list(values = subarray(reg_births, time == "2011-2015"), col = "red", lwd= 2))

dplot( ~ age | region, data = dea.chain,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Deaths estimation 2006-2014",
       overlay = list(values = subarray(reg_deaths, time == "2011-2015"), col = "red", lwd= 2))

dplot( ~ age | region, data = imm.chain,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Immigration estimation 2006-2014",
       overlay = list(values = subarray(arrivals, time == "2011-2015"), col = "red", lwd= 2))

dplot( ~ age | region, data = emi.chain,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Emigration estimation 2006-2014",
       overlay = list(values = subarray(departures, time == "2011-2015"), col = "red", lwd= 2))

dplot( ~ age | region, data = iin.chain ,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal immigration estimation 2006-2014",
       overlay = list(values =subarray(intarr, time == "2011-2015"), col = "red", lwd= 2))

dplot( ~ age | region, data = ein.chain,
       subarray=time=="2011-2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal emigration estimation 2006-2014",
       overlay = list(values = subarray(intdep, time == "2011-2015"), col = "red", lwd= 2))


MCMC <- fetchMCMC(onlyreg)
names(MCMC)
library(coda)
cbind(
  dimnames(MCMC$account.population[[1]])[[2]],
  dimnames(MCMC$account.population[[2]])[[2]],
  dimnames(MCMC$account.population[[3]])[[2]],
  dimnames(MCMC$account.population[[4]])[[2]]
)


plot(MCMC$account.population[[1]])
dim(MCMC$systemModels.births.likelihood.rate[[1]])
dimnames(MCMC$systemModels.births.likelihood.rate[[1]])
plot(MCMC$systemModels.births.likelihood.rate[[1]])
plot(MCMC$systemModels.births.prior.mean[[1]])

gelman.diag(MCMC$systemModels.population.likelihood.count)
