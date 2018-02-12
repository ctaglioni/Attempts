#-------------------
# Only age
#-------------------

library(methods)
#library(nzreg)
library(dplyr)
library(docopt)
library(demest)
library(demItaly)
library(abind)
library(beepr)

dimnames(italy.popn.reg)$time <- 2001:2016
pop <- italy.popn.reg[,,-c(5,6),]
bir <- italy.births.reg[-c(5,6),,]
dea <- italy.deaths.reg[,,-c(5,6),]
ei <- italy.ext.in[,,-c(5,6),]
ee <- italy.ext.out[,,-c(5,6),]

census_year_erp<- Counts(pop, dimscales = c(time="Points")) %>%
  collapseDimension(margin = c("age", "time"))
# Need 2005 as the scale is Point

# as 5years age classes only take years 2006-2010 
# (avoid error:'population' does not have regular age-time plan : age step [5] does not equal time step [1])
population<- census_year_erp %>%
  subarray( time %in% c("2005", "2010", "2015"))
population <- Counts(population, dimscales = c(time="Points"))

reg_births<- Counts(bir, dimscales = c(time="Intervals")) %>%
  subarray(time != "2016") %>%
  reallocateToEndAges() %>%
  collapseDimension(margin = c("age", "time")) %>%
  collapseIntervals( dimension= "time", width = 5) %>%
  collapseIntervals(dimension = "age", width = 5)
# plot(reg_births)

reg_deaths <- Counts(dea, dimscales = c(time="Intervals")) %>%
  subarray(time != "2016") %>%
  collapseDimension(margin = c("age", "time")) %>%
  collapseIntervals(dimension = "time", width = 5)
# plot(reg_deaths)

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

arr <- MigrAgeGroupDisag(ei, r=20)
dim(arr)
dimnames(arr) <-list(age=dimnames(pop)[[1]],
                     sex=dimnames(pop)[[2]],
                     region = dimnames(pop)[[3]],
                     time = c(2006:2016))

arrivals <- Counts(arr , dimscales = c(time="Intervals")) %>%
  subarray(time != "2016") %>%
  #expandIntervals(dimension = "age", width = 5, old = dimnames(italy.ext.in)[[1]])
  collapseDimension(margin = c("age", "time")) %>%
  collapseIntervals(dimension = "time", width = 5)
# plot(arrivals)


dep <-MigrAgeGroupDisag(ee, r=20)
dim(dep)
dimnames(dep) <-list(age=dimnames(pop)[[1]],
                     sex=dimnames(pop)[[2]],
                     region = dimnames(pop)[[3]],
                     time = c(2006:2016))

departures <- Counts(dep, dimscales = c(time="Intervals")) %>%
  subarray(time != "2016") %>%
  collapseDimension(margin = c("age", "time"))%>%
  collapseIntervals(dimension = "time", width = 5)
# plot(departures)

account<- Movements(population = population,
                    births = reg_births,
                    entries = list(external_in = arrivals),
                    exits = list(deaths = reg_deaths,
                                 external_out = departures)) %>%
  makeConsistent()

# summary(account)
systemModels<- list(Model(population ~ Poisson(mean ~ age + time + age*time, useExpose=FALSE),
                          time ~ DLM(),
                          age ~ Exch(),
                          jump = 0.002),
                    Model(births ~ Poisson(mean ~ 1),
                          jump = 0.005),
                    Model(deaths ~ Poisson(mean ~ 1),
                          jump = 0.01),
                    Model(external_in ~ Poisson(mean ~ 1),
                          jump = 0.01),
                    Model(external_out ~ Poisson(mean ~ 1),
                          jump = 0.03))

datasets <- list(census_year_erp = population,
                 reg_births = reg_births,
                 reg_deaths = reg_deaths,
                 arrivals = arrivals,
                 departures = departures)

sd<- census_year_erp %>%
  collapseDimension(margin = c("age", "time")) %>%
  as("Values") 
sd<- sd * 0.0025
mean<- sd / sd

sdd <- reg_deaths %>%
  collapseDimension(margin = c("age", "time")) %>%
  as("Values") 
sdd<- sdd * 0.0025
meand<- sdd / sdd

dataModels <- list(Model(census_year_erp ~ NormalFixed(mean = mean, sd = sd),
                         series = "population"),
                   Model(reg_births ~ PoissonBinomial(prob = 0.95),
                         series = "births"),
                   Model(reg_deaths ~ NormalFixed(mean = meand, sd = sdd),
                         series = "deaths"),
                   Model(arrivals ~ Poisson(mean ~ 1),
                         series = "external_in",
                         jump = 0.01),
                   Model(departures ~ Poisson(mean ~ 1),
                         series = "external_out",
                         jump = 0.01))


onlyage <- "C:\\0_PhD\\Thesis\\Thesis_R\\AT9.est"

n_sim <- 50000
n_burnin <- 50000
n_chain <- 4
n_thin <- 400

beep(sound=2, estimateAccount(account = account,
                systemModels = systemModels,
                datasets = datasets,
                dataModels = dataModels,
                filename = onlyage,
                nBurnin = n_burnin,
                nSim = n_sim,
                nChain = n_chain,
                nThin = n_thin,
                useC = TRUE))

fetchSummary(onlyage)
showModel(onlyage)
listContents(onlyage)

pop.chain<- fetch(onlyage, where=c("account", "population"))
bir.chain <- fetch(onlyage, where=c("account", "births"))
dea.chain <- fetch(onlyage, where=c("account", "deaths"))
imm.chain <- fetch(onlyage, where=c("account", "external_in"))
emi.chain <- fetch(onlyage, where=c("account", "external_out"))


firstrow <- 0.5*(n_sim/n_thin*n_chain)
lastrow <- dim(pop.chain)[3]

dim(pop.chain[,,firstrow:lastrow])
popest<- apply(pop.chain[,,firstrow:lastrow], c(1,2), mean)
birest2<- apply(bir.chain[,,,firstrow:lastrow], c(1,2,3), mean)
birest<- apply(birest2, c(1,2),sum)
deaest2 <- apply(dea.chain[,,,firstrow:lastrow], c(1,2,3), mean)
deaest<- apply(deaest2, c(1,2),sum)
immest2 <- apply(imm.chain[,,,firstrow:lastrow], c(1,2,3), mean)
immest<- apply(immest2, c(1,2),sum)
emiest2 <- apply(emi.chain[,,,firstrow:lastrow], c(1,2,3), mean)
emiest<- apply(emiest2, c(1,2),sum)

# Estimate - truth
dplot( ~ age | time, data = popest - population,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Population differences 2006-2015")

plot(popest[,1] - population[,1], type="l", xaxt="n",
     xlab = "Age", ylab="Estimate - Data", 
     main = "Population: estimate vs data")
axis(1, at = 1:19 , labels = dimnames(popest)[[1]])
lines(popest[,2] - population[,2], col = 2)
lines(popest[,3] - population[,3], col = 3)

dplot( ~ age | time, data = birest - reg_births,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births differences 2006-2015")
dev.off()
plot(birest[,1] - reg_births[,1], type="l", xaxt="n",
     xlab = "Age", ylab="Estimate - Data", 
     main = "Births: estimate vs data",
     ylim = c(min(birest - reg_births), max(birest - reg_births)))
axis(1, at = 1:7 , labels = dimnames(birest)[[1]])
lines(birest[,2] - reg_births[,2], col = 2)

dplot( ~ age | time, data = deaest - reg_deaths,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Deaths differences 2006-2015")
plot(deaest[,1] - reg_deaths[,1], type="l", xaxt="n",
     xlab = "Age", ylab="Estimate - Data", 
     main = "Deaths: estimate vs data",
     ylim = c(min(deaest - reg_deaths), max(deaest - reg_deaths)))
axis(1, at = 1:19 , labels = dimnames(deaest)[[1]])
lines(deaest[,2] - reg_deaths[,2], col = 2)


dplot( ~ age | time, data = immest - arrivals,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Immigration differences 2006-2015")
plot(immest[,1] - arrivals[,1], type="l", xaxt="n",
     xlab = "Age", ylab="Estimate - Data", 
     main = "Immigration: estimate vs data",
     ylim = c(min(immest - arrivals), max(immest - arrivals)))
axis(1, at = 1:19 , labels = dimnames(arrivals)[[1]])
lines(immest[,2] - arrivals[,2], col = 2)

dplot( ~ age | time, data = emiest - departures,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Emigration differences 2006-2015")
plot(emiest[,1] - departures[,1], type="l", xaxt="n",
     xlab = "Age", ylab="Estimate - Data", 
     main = "Emigration: estimate vs data",
     ylim = c(min(emiest - departures), max(emiest - departures)))
axis(1, at = 1:19 , labels = dimnames(departures)[[1]])
lines(emiest[,2] - departures[,2], col = 2)


#----------------------

dplot( ~ age | time, data = pop.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Population estimation 2006-2014",
       overlay = list(values = population, col = "red", lwd= 2))


dplot( ~ age | time, data = bir.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation 2006-2014",
       overlay = list(values = reg_births, col = "red", lwd= 2))

dplot( ~ age | time, data = dea.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Deaths estimation 2006-2014",
       overlay = list(values = reg_deaths, col = "red", lwd= 2))

dplot( ~ age | time, data = imm.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Immigration estimation 2006-2014",
       overlay = list(values = arrivals, col = "red", lwd= 2))

dplot( ~ age | time, data = emi.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Emigration estimation 2006-2014",
       overlay = list(values = departures, col = "red", lwd= 2))


MCMC <- fetchMCMC(onlyage)
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
plot(MCMC$systemModels.deaths.hyper.time.scaleError, sub=names(MCMC)[22])
plot(MCMC$systemModels.external_out.likelihood.rate, sub=names(MCMC)[23])
plot(MCMC$systemModels.external_out.prior.mean, sub=names(MCMC)[24]) # 2014! But seems to converge at the end
plot(MCMC$systemModels.external_out.hyper.time.scaleLevel, sub=names(MCMC)[25]) # prior ok
plot(MCMC$systemModels.external_out.hyper.time.scaleError, sub=names(MCMC)[26])
plot(MCMC$dataModels.arrivals.likelihood.rate, sub=names(MCMC)[27]) # as system.model
plot(MCMC$dataModels.arrivals.prior.mean, sub=names(MCMC)[28])
plot(MCMC$dataModels.arrivals.prior.sd, sub=names(MCMC)[29]) #...convergence?
plot(MCMC$dataModels.departures.likelihood.rate, sub=names(MCMC)[30]) # most problematic chain
plot(MCMC$dataModels.departures.prior.mean, sub=names(MCMC)[31]) # prior ok
plot(MCMC$dataModels.departures.prior.sd, sub=names(MCMC)[32])
