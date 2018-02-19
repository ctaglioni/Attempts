#-------------------
# Only region
#-------------------

library(methods)
#library(nzreg)
library(dplyr)
library(docopt)
library(demest)
library(demItaly)
library(beepr)

dimnames(italy.popn.reg)$time <- 2001:2016
census_year_erp<- Counts(italy.popn.reg, dimscales = c(time="Points")) %>%
  collapseDimension(margin = c("time", "region"))%>%
  subarray(time != c("2016")) %>%
  subarray(time >"2004")

population<- census_year_erp

reg_births<- Counts(italy.births.reg, dimscales = c(time="Intervals")) %>%
  subarray(time !="2016") %>%
  collapseDimension(margin = c("time", "region"))
# plot(reg_births)

reg_deaths <- Counts(italy.deaths.reg, dimscales = c(time="Intervals")) %>%
  collapseDimension(margin = c("time", "region"))
# plot(reg_deaths)

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

systemModels<- list(Model(population ~ Poisson(mean ~ time * region, useExpose=FALSE),
                          time : region ~ Exch(),
                          jump = 0.002),
                    Model(births ~ Poisson(mean ~ time * region),
                          time : region ~ Exch(),
                          jump = 0.015),
                    Model(deaths ~ Poisson(mean ~ time * region),
                          time : region ~ Exch(),
                          jump = 0.015),
                    Model(external_in ~ Poisson(mean ~ time * region),
                          time : region ~ Exch(),
                          jump = 0.03),
                    Model(external_out ~ Poisson(mean ~ time * region),
                          time : region ~ Exch(),
                          jump = 0.03),
                    Model(internal_in ~ Poisson(mean ~ time * region),
                          time : region ~ Exch(),
                          jump = 0.03),
                    Model(internal_out ~ Poisson(mean ~ time * region),
                          time : region ~ Exch(),
                          jump = 0.03))

# the population size of the region is not directly linked to the number of migrants
# i.e. only exchangeable prior, nor covariate

cy1 <- census_year_erp %>%
  subarray(time != "2011")

cy2 <- census_year_erp%>%
  subarray(time == "2011")
cy2 <- as.array(cy2)
dimnames(cy2)<- list(region = dimnames(cy1)[[2]])
cy2 <- Counts(cy2)
cy2011 <- addDimension(cy2, name = "time", labels = "2011", dimscale = c(time ="Points"))


datasets <- list(census_year_erp1 = cy1,
                 census_year_erp2 = cy2011,
                 reg_births = reg_births,
                 reg_deaths = reg_deaths,
                 arrivals = arrivals,
                 departures = departures,
                 intarr = intarr,
                 intdep = intdep)



dataModels <- list(Model(census_year_erp1 ~ Poisson(mean ~ region),
                         series = "population",
                         jump = 0.002),
                   Model(census_year_erp2 ~ PoissonBinomial(prob = 0.95),
                         series = "population"),
                   Model(reg_births ~ PoissonBinomial(prob = 0.95),
                         series = "births"),
                   Model(reg_deaths ~ PoissonBinomial(prob = 0.90),
                         series = "deaths"),
                   Model(arrivals ~ Poisson(mean ~ region),
                         region ~ Exch(covariates = Covariates(mean ~ mean, data = mean_pop)),
                         series = "external_in",
                         jump = 0.03),
                   Model(departures ~ Poisson(mean ~ region),
                         region ~ Exch(covariates = Covariates(mean ~ mean, data = mean_pop)),
                         series = "external_out",
                         jump = 0.03),
                   Model(intarr ~ Poisson(mean ~ region),
                         region ~ Exch(covariates = Covariates(mean ~ mean, data = mean_pop)),
                         series = "internal_in",
                         jump = 0.03),
                   Model(intdep ~ Poisson(mean ~ region),
                         region ~ Exch(covariates = Covariates(mean ~ mean, data = mean_pop)),
                         series = "internal_out",
                         jump = 0.03))


onlyreg <- "C:/0_PhD/Thesis/Thesis_R/RT5.est"

n_sim <- 50000
n_burnin <- 50000
n_chain <- 3
n_thin <- 250

beep(sound=2, estimateAccount(account = account,
                systemModels = systemModels,
                datasets = datasets,
                dataModels = dataModels,
                filename = onlyreg,
                nBurnin = n_burnin,
                nSim = n_sim,
                nChain = n_chain,
                nThin = n_thin,
                useC = TRUE))

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

par(mfrow=c(2,4))
plot(popest[11,] - population[11,], main = "Population: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:22 , labels = dimnames(popest)[[2]], cex.axis=0.4, las = 3)
abline(h=0)


plot(birest[10,] - reg_births[10,], main = "Births: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:22 , labels = dimnames(birest)[[2]], cex.axis=0.4, las = 3)
abline(h=0)

plot(deaest[10,] - reg_deaths[10,], main = "Deaths: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:22 , labels = dimnames(deaest)[[2]], cex.axis=0.4, las = 3)
abline(h=0)

plot(immest[10,] - arrivals[10,], main = "Ext immigration: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:22 , labels = dimnames(immest)[[2]], cex.axis=0.4, las = 3)
abline(h=0)

plot(emiest[10,] - departures[10,], main = "Ext emigration: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:22 , labels = dimnames(emiest)[[2]], cex.axis=0.4, las = 3)
abline(h=0)

plot(iinest[10,] - intarr[10,], main = "Int immigration: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:22 , labels = dimnames(iinest)[[2]], cex.axis=0.4, las = 3)
abline(h=0)

plot(einest[10,] - intdep[10,], main = "Int emigration: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:22 , labels = dimnames(einest)[[2]], cex.axis=0.4, las = 3)
abline(h=0)

par(mfrow=c(1,1))
plot(popest[2:11,1] - population[2:11,1], main = "Population: estimation vs data",
     type = "l", ylim=c(min(popest[2:11,] - population[2:11,]),
                        max(popest[2:11,] - population[2:11,])),
     xlab = "Year", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:10 , labels = dimnames(popest)[[1]][2:11])
abline(h=0)
for (i in 2:22){
  lines(popest[2:11,i] - population[2:11,i], col=i)
}


plot(immest[2:10,1] - arrivals [2:10,1], main = "External immigration",
     type = "l", ylim=c(min(immest[2:10,] - arrivals[2:10,]),
                        max(immest[2:10,] - arrivals[2:10,])),
     xlab = "Year", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:10 , labels = dimnames(popest)[[1]][2:11])
abline(h=0)
for (i in 2:22){
  lines(immest[2:10,i] - arrivals[2:10,i], col=i)
}
legend("topright", dimnames(immest)[[2]], lty=1, col = 2:22, cex = 0.5)


plot(emiest[2:10,1] - departures [2:10,1], main = "External emigration",
     type = "l", ylim=c(min(emiest[2:10,] - departures[2:10,]),
                        max(emiest[2:10,] - departures[2:10,])),
     xlab = "Year", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:10 , labels = dimnames(popest)[[1]][2:11])
abline(h=0)
for (i in 2:22){
  lines(emiest[2:10,i] - departures[2:10,i], col=i)
}
legend("topright", dimnames(popest)[[2]], lty=1, col = 2:22, cex = 0.5)

plot(iinest[2:10,1] - intarr [2:10,1], main = "Internal immigration",
     type = "l", ylim=c(min(iinest[2:10,] - intarr[2:10,]),
                        max(iinest[2:10,] - intarr[2:10,])),
     xlab = "Year", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:10 , labels = dimnames(popest)[[1]][2:11])
abline(h=0)
for (i in 2:22){
  lines(iinest[2:10,i] - intarr[2:10,i], col=i)
}
legend("topright", dimnames(popest)[[2]], lty=1, col = 2:22, cex = 0.5)

plot(apply(iinest[2:10,] - intarr [2:10,], 2 , mean), apply(population[-1,], 2, mean))


plot(einest[2:10,1] - intdep [2:10,1], main = "Internal emigration",
     type = "l", ylim=c(min(einest[2:10,] - intdep[2:10,]),
                        max(einest[2:10,] - intdep[2:10,])),
     xlab = "Year", ylab = "Diff", xaxt = "n", pch = 18)
axis(1,at = 1:10 , labels = dimnames(popest)[[1]][2:11])
abline(h=0)
for (i in 2:22){
  lines(einest[2:10,i] - intdep[2:10,i], col=i)
}
legend("topright", dimnames(popest)[[2]], lty=1, col = 2:22, cex = 0.5)

plot(apply(einest[2:10,] - intdep [2:10,], 2 , mean), apply(population[-1,], 2, mean))



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
       overlay = list(values = arrivals, col = "red", lwd= 2))

dplot( ~ time | region, data = ein.chain,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Internal emigration estimation 2006-2014",
       overlay = list(values = departures, col = "red", lwd= 2))


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

