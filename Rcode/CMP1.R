# Dispersion(mean = 0.1, scale = HalfT(scale = 0.05))
library(dplyr)
library(demest)
library(demItaly)
library(demCMP)
library(beepr)


#---------------------------------------------------
# Simulated underdispersed data coverage regression
#---------------------------------------------------
age <- round(rnorm(5,0,0.1),3) # 0
A <- array(age, dim=c(5, 2, 10, 3))
sex <- cbind(rep(0,5), rep(-0.1,5)) # na
S <- array(sex,dim=c(5, 2, 10, 3))
region <- round(rnorm(10,-1,0.1), 3) # -1
R <- array(NA, dim=c(5, 2, 10,3))
for(i in 1:length(region)){
  R[,,i,] <- matrix(region[i], 5,2)
}
time <- round(rnorm(3,0.5,0.1), 3) #0.5
Y <- array(NA, dim=c(5, 2, 10,3))
for(i in 1:length(time)){
  Y[,,,i] <- array(time[i], c(5,2,10))
}

mu <- R+S+A+Y

sigma <- 0.1
eta <- 0.5
tau <- 0.05
expo <- replicate(300, rcmp1(mu = 1000, nu = 1 , max = 10000))
gamma<- exp(rnorm(300, mean = mu, sd = sigma))
nu <- exp(rnorm(300, mean = eta, sd = tau))
mu.cmp<- expo*gamma
y <- c()
for(i in 1:length(expo)){
  y[i] <- rcmp1(mu=mu.cmp[i], nu=nu[i], max=1000)
}

expodem <- array(expo, dim=c(5, 2, 10, 3))
dimnames(expodem) <- list(age = 1:5,
                          sex = c("f","m"),
                          region = letters[1:10],
                          time = 2014:2016)

ydem <- array(y, dim=c(5, 2, 10, 3))
dimnames(ydem) <- list(age = 1:5,
                       sex = c("f","m"),
                       region = letters[1:10],
                       time = 2014:2016)

yy <- Counts(ydem, dimscales = c(time ="Intervals"))
exx <- Counts(expodem, dimscales = c(time ="Intervals"))

filename <- tempfile()

beep(sound=2, estimateModel(Model(yy ~ CMP(mean ~ age +sex + region + time,
                                           dispersion = Dispersion(mean = Norm(mean = 0, sd = 1),
                                                                   scale = HalfT(scale = 0.5))),
                                  age ~ Exch(),
                                  region ~ Exch(),
                                  time ~ Exch(),
                                  priorSD = HalfT(scale = 0.5),
                                  jump = 0.03),
                            y = yy,
                            exposure = exx,
                            filename = filename,
                            nBurnin = 100000,
                            nSim = 100000,
                            nChain = 3,
                            nThin = 100))

beep(sound=2, estimateModel(Model(yy ~ CMP(mean ~ 1,
                                           dispersion = Dispersion(mean = Norm(mean = 0, sd = 1),
                                                                   scale = HalfT(scale = 0.5))),
                                  jump = 0.03),
                            y = yy,
                            exposure = exx,
                            filename = filename,
                            nBurnin = 100000,
                            nSim = 100000,
                            nChain = 3,
                            nThin = 100))


beep(sound=2,continueEstimation(filename,nBurnin = 300000,
                                nSim = 300000,
                                nThin = 1000))

fetchSummary(filename)
showModel(filename)
listContents(filename)

rate.lik <- fetch(filename, where=c("model", "likelihood", "rate"))
disp.lik <- fetch(filename, where=c("model", "likelihood", "dispersion"))
disp.prior <- fetch(filename, where=c("model", "prior", "dispersion", "mean"))


rateest<- apply(rate.lik, c(1,2,3,4), mean)
dispest <- apply(disp.lik, c(1,2,3,4), mean)
par(mfrow=c(1,1))
plot(rateest[,1,1,1]*expodem[,1,1,1] - yy[,1,1,1], main = "estimation vs data 2015",
     xlab = "Age", ylab = "Diff", pch = 18)


dplot( ~ age | region, data = rate.lik[,1,,,], subarray = time == "2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Rate estimation")

dplot( ~ age | region, data = disp.lik[,1,,,], subarray = time == "2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Dispersion estimation")

plot(rateest - array(gamma, dim=c(5, 2, 10, 3)))
abline(h=0)

plot(dispest - array(nu, dim=c(5, 2, 10, 3)))
abline(h=0)

MCMC <-fetchMCMC(filename)

plot(MCMC$model.likelihood.rate, sub= "Likelihood rate")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.prior.rate.mean, sub= "Prior rate mean" )
plot(MCMC$model.prior.rate.sd, sub= "Prior rate sd" )







#---------------------------------------------------
# Simulated underdispersed data coverage more random
#---------------------------------------------------

mu <- -0.1
sigma <- 0.1
eta <- -0.2
tau <- 0.05
expo <- replicate(2000, rcmp1(mu = 1000, nu = 1 , max = 10000))
gamma<- exp(rnorm(2000, mean = mu, sd = sigma))
nu <- exp(rnorm(2000, mean = eta, sd = tau))
mu.cmp<- expo*gamma
y <- c()
for(i in 1:length(expo)){
  y[i] <- rcmp1(mu=mu.cmp[i], nu=nu[i], max=1000)
}

expodem <- array(expo, dim=c(5, 2, 10, 2))
dimnames(expodem) <- list(age = 1:5,
                          sex = c("f","m"),
                          region = letters[1:10],
                          time = 2015:2016)

ydem <- array(y, dim=c(5, 2, 10, 2))
dimnames(ydem) <- list(age = 1:5,
                       sex = c("f","m"),
                       region = letters[1:10],
                       time = 2015:2016)

yy <- Counts(ydem, dimscales = c(time ="Intervals"))
exx <- Counts(expodem, dimscales = c(time ="Intervals"))

filename <- tempfile()

beep(sound=2, estimateModel(Model(yy ~ CMP(mean ~ 1,
                             dispersion = Dispersion(mean = Norm(mean = -0.25, sd = 0.07),
                                                     scale = HalfT(scale = 0.3))),
                    jump = 0.03),
              y = yy,
              exposure = exx,
              filename = filename,
              nBurnin = 500000,
              nSim = 500000,
              nChain = 3,
              nThin = 1000))

 beep(sound=2,continueEstimation(filename,nBurnin = 300000,
                   nSim = 300000,
                   nThin = 1000))

fetchSummary(filename)
showModel(filename)
listContents(filename)

rate.lik <- fetch(filename, where=c("model", "likelihood", "rate"))
disp.lik <- fetch(filename, where=c("model", "likelihood", "dispersion"))
disp.prior <- fetch(filename, where=c("model", "prior", "dispersion", "mean"))


rateest<- apply(rate.lik, c(1,2,3,4), mean)
dispest <- apply(disp.lik, c(1,2,3,4), mean)
par(mfrow=c(1,1))
plot(rateest[,1,1,1]*expodem[,1,1,1] - yy[,1,1,1], main = "estimation vs data 2015",
     xlab = "Age", ylab = "Diff", pch = 18)


dplot( ~ age | region, data = rate.lik[,1,,,], subarray = time == "2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Rate estimation")

dplot( ~ age | region, data = disp.lik[,1,,,], subarray = time == "2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Dispersion estimation")

plot(rateest - array(gamma, dim=c(5, 2, 10, 2)))
abline(h=0)

plot(dispest - array(nu, dim=c(5, 2, 10, 2)))
abline(h=0)

MCMC <-fetchMCMC(filename)

plot(MCMC$model.likelihood.rate, sub= "Likelihood rate")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.prior.rate.mean, sub= "Prior rate mean" )
plot(MCMC$model.prior.rate.sd, sub= "Prior rate sd" )




#------------------------------------------
# Simulated underdispersed data coverage ~1
#------------------------------------------

mu <- 0
sigma <- 0.01
eta <- 1
tau <- 0.05
expo <- replicate(200, rcmp1(mu = 1000, nu = 1 , max = 10000))
gamma <- exp(rnorm(200, mean = mu, sd = sigma))
nu <- exp(rnorm(200, mean = eta, sd = tau))
mu.cmp<- expo*gamma
y <- c()
for(i in 1:length(expo)){
  y[i] <- rcmp1(mu=mu.cmp[i], nu=nu[i], max=1000)
}

expodem <- array(expo, dim=c(5, 2, 10, 2))
dimnames(expodem) <- list(age = 1:5,
                          sex = c("f","m"),
                          region = letters[1:10],
                          time = 2015:2016)

ydem <- array(y, dim=c(5, 2, 10, 2))
dimnames(ydem) <- list(age = 1:5,
                       sex = c("f","m"),
                       region = letters[1:10],
                       time = 2015:2016)

yy <- Counts(ydem, dimscales = c(time ="Intervals"))
exx <- Counts(expodem, dimscales = c(time ="Intervals"))

filename <- tempfile()

estimateModel(Model(yy ~ CMP(mean ~ 1,
                             dispersion = Dispersion(mean = Norm(mean = 0, sd = 1),
                                                     scale = HalfT(scale = 0.1))),
                    jump = 0.01),
              y = yy,
              exposure = exx,
              filename = filename,
              nBurnin = 500000,
              nSim = 500000,
              nChain = 3,
              nThin = 1000)

# continueEstimation(filename,nBurnin = 300000,
#                    nSim = 300000,
#                    nThin = 1000)

fetchSummary(filename)
showModel(filename)
listContents(filename)

rate.lik <- fetch(filename, where=c("model", "likelihood", "rate"))
disp.lik <- fetch(filename, where=c("model", "likelihood", "dispersion"))
disp.prior <- fetch(filename, where=c("model", "prior", "dispersion", "mean"))


rateest<- apply(rate.lik[,,,,800:1500], c(1,2,3,4), mean)
dispest <- apply(disp.lik[,,,,800:1500], c(1,2,3,4), mean)

plot(rateest[,1,1,1]*expodem[,1,1,1] - yy[,1,1,1], main = "estimation vs data 2015",
     xlab = "Age", ylab = "Diff", pch = 18)


dplot( ~ age | region, data = rate.lik[,1,,,], subarray = time == "2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Rate estimation")

dplot( ~ age | region, data = disp.lik[,1,,,], subarray = time == "2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Dispersion estimation")

plot(rateest - array(gamma, dim=c(5, 2, 10, 2)))
abline(h=0)

plot(dispest - array(nu, dim=c(5, 2, 10, 2)))
abline(h=0)

MCMC <-fetchMCMC(filename)

plot(MCMC$model.likelihood.rate, sub= "Likelihood rate")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.prior.rate.mean, sub= "Prior rate mean" )
plot(MCMC$model.prior.rate.sd, sub= "Prior rate sd" )



#--------------------------------
# >Simulated underdispersed data
#--------------------------------

mu <- -0.2
sigma <- 0.01
eta <- 1
tau <- 0.05
expo <- replicate(200, rcmp1(mu = 1000, nu = 0.6 , max = 1000))
gamma <- exp(rnorm(200, mean = mu, sd = sqrt(sigma)))
nu<- exp(rnorm(200, mean = eta, sd = sqrt(tau)))
mu.cmp<- expo*gamma
y <- c()
for(i in 1:length(expo)){
  y[i] <- rcmp1(mu=mu.cmp[i], nu=nu[i], max=1000)
}
y

expodem <- array(expo, dim=c(5, 2, 10, 2))
dimnames(expodem) <- list(age = 1:5,
                          sex = c("f","m"),
                          region = letters[1:10],
                          time = 2015:2016)

ydem <- array(y, dim=c(5, 2, 10, 2))
dimnames(ydem) <- list(age = 1:5,
                          sex = c("f","m"),
                          region = letters[1:10],
                          time = 2015:2016)

yy <- Counts(ydem, dimscales = c(time ="Intervals"))
exx <- Counts(expodem, dimscales = c(time ="Intervals"))

filename <- tempfile()

estimateModel(Model(yy ~ CMP(mean ~ age + region + time,
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 0.1),
                                                    scale = HalfT(scale = 0.1))),
                    jump = 0.02),
              y = yy,
              exposure = exx,
              filename = filename,
              nBurnin = 50000,
              nSim = 50000,
              nChain = 3,
              nThin = 200)

continueEstimation(filename,nBurnin = 30000,
                   nSim = 30000,
                   nThin = 750)

fetchSummary(filename)
showModel(filename)
listContents(filename)

bir.chain <- fetch(filename, where=c("model", "likelihood", "rate"))
disp.lik <- fetch(filename, where=c("model", "likelihood", "dispersion"))
disp.prior <- fetch(filename, where=c("model", "prior", "dispersion", "mean"))


birest<- apply(bir.chain, c(1,2,3,4), mean)
dispest <- apply(disp.lik, c(1,2,3,4), mean)

plot(birest[,,1,1]*expodem[,,1,1] - yy[,,1,1], main = "Births: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", pch = 18)

dplot( ~ age | region, data = bir.chain[,1,,,], subarray = time == "2015",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation")

plot(birest - array(gamma, dim=c(5, 2, 10, 2)))
abline(h=0)

plot(dispest - array(nu, dim=c(5, 2, 10, 2)))
abline(h=0)

data <- Counts(birest[,,11]*expose[,,15])
dplot( ~ age | region, data = data,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation",
       overlay = list(values = Counts(italy.births.reg[,,11]), col = "red", lwd= 2))


MCMC <-fetchMCMC(filename)

plot(MCMC$model.likelihood.rate, sub= "Likelihood rate")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.prior.rate.mean, sub= "Prior rate mean" )
plot(MCMC$model.prior.rate.sd, sub= "Prior rate sd" )



#------------------------------
# Deaths without exposure
#------------------------------

y <- italy.deaths.reg %>%
  Counts(dimscales = c(time = "Intervals"))

filename <- tempfile()

estimateModel(Model(y ~ CMP(mean ~ age + time +region,
                            dispersion = Dispersion(mean = Norm(mean = 0
                                                                , sd = 1),
                                                    scale = HalfT(scale = 0.1)),
                            useExpose = FALSE),
                    jump = 0.03),
              y = y,
              filename = filename,
              nBurnin = 50000,
              nSim = 50000,
              nChain = 3,
              nThin = 250)

continueEstimation(filename,nBurnin = 60000,
                   nSim = 60000,
                   nThin = 750)

fetchSummary(filename)
showModel(filename)
listContents(filename)

dea.chain <- fetch(filename, where=c("model", "likelihood", "count"))
disp.lik <- fetch(filename, where=c("model", "likelihood", "dispersion"))
disp.prior <- fetch(filename, where=c("model", "prior", "dispersion", "mean"))


deaest<- apply(dea.chain, c(1,2,3), mean)
dispest<- apply(disp.lik, c(1,2,3), mean)
par(mfrow=c(1,1))
plot(deaest[,7,11] - italy.deaths.reg[,7,11], main = "deaths: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", pch = 18)

dplot( ~ age | region, data = dea.chain, subarray = time == "2012",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "deaths estimation")

dplot( ~ age | region, data = disp.lik, subarray = time == "2012",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "deaths estimation")


dplot( ~ age | region, data = Counts(deaest, dimscales = c(time = "Intervals")),
       subarray = time == "2012",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "deaths estimation",
       overlay = list(values = Counts(italy.deaths.reg[,,11]), col = "red", lwd= 2))


MCMC <-fetchMCMC(filename)

plot(MCMC$model.likelihood.count, sub= "Likelihood counts")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.hyper.age.scaleLevel, sub= "Hyper age scale level" )
plot(MCMC$model.prior.count.mean, sub= "Prior count mean" )
plot(MCMC$model.prior.count.sd, sub= "Prior count sd" )
plot(MCMC$model.hyper.age.scaleTrend, sub= "Hyper age scale trend")
plot(MCMC$model.hyper.age.scaleError, sub= "Hyper age scale error" )
plot(MCMC$model.hyper.year.damp, sub= "Hyper year damp" )
plot(MCMC$model.hyper.year.scaleError, sub= "Hyper year scale error" )

#------------------------------
# Deaths with exposure
#------------------------------
y <- italy.deaths.reg %>%
  Counts(dimscales = c(time = "Intervals"))

expose <- italy.popn.reg %>%
  Counts(dimscales = c(time = "Intervals"))

filename <- tempfile()

estimateModel(Model(y ~ CMP(mean ~ age + region + time,
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 1),
                                                    scale = HalfT(scale = 0.1))),
                    age ~ DLM(damp = NULL),
                    jump = 0.02),
              y = y,
              exposure = expose,
              filename = filename,
              nBurnin = 300000,
              nSim = 300000,
              nChain = 3,
              nThin = 1000)

continueEstimation(filename,nBurnin = 30000,
                   nSim = 30000,
                   nThin = 750)

fetchSummary(filename)
showModel(filename)
listContents(filename)

dea.chain <- fetch(filename, where=c("model", "likelihood", "rate"))
disp.lik <- fetch(filename, where=c("model", "likelihood", "dispersion"))
disp.prior <- fetch(filename, where=c("model", "prior", "dispersion", "mean"))


deaest<- apply(dea.chain, c(1,2,3), mean)
dispest <- apply(disp.lik, c(1,2,3), mean)

expose <- aperm(expose, c(2,1,3))
plot(deaest[,,11]*expose[,,15] - italy.deaths.reg[,,11], main = "deaths: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", pch = 18)

dplot( ~ age | region, data = dea.chain, subarray = time == "2012",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "deaths estimation")

data <- Counts(deaest[,,11]*expose[,,15])
dplot( ~ age | region, data = data,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "deaths estimation",
       overlay = list(values = Counts(italy.deaths.reg[,,11]), col = "red", lwd= 2))


MCMC <-fetchMCMC(filename)

plot(MCMC$model.likelihood.rate, sub= "Likelihood rate")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.hyper.age.scaleLevel, sub= "Hyper age scale level" )
plot(MCMC$model.prior.rate.mean, sub= "Prior rate mean" )
plot(MCMC$model.prior.rate.sd, sub= "Prior rate sd" )
plot(MCMC$model.hyper.age.scaleTrend, sub= "Hyper age scale trend")
plot(MCMC$model.hyper.age.scaleError, sub= "Hyper age scale error" )





#------------------------------
# Births without exposure
#------------------------------

y <- italy.births.reg %>%
  Counts(dimscales = c(time = "Intervals"))

filename <- tempfile()

estimateModel(Model(y ~ CMP(mean ~ age + time +region,
                            dispersion = Dispersion(mean = Norm(mean = 0
                                                                , sd = 1),
                                                    scale = HalfT(scale = 0.1)),
                            useExpose = FALSE),
                    jump = 0.02),
              y = y,
              filename = filename,
              nBurnin = 300000,
              nSim = 300000,
              nChain = 3,
              nThin = 1000)

continueEstimation(filename,nBurnin = 60000,
                   nSim = 60000,
                   nThin = 750)

fetchSummary(filename)
showModel(filename)
listContents(filename)

bir.chain <- fetch(filename, where=c("model", "likelihood", "count"))
disp.lik <- fetch(filename, where=c("model", "likelihood", "dispersion"))
disp.prior <- fetch(filename, where=c("model", "prior", "dispersion", "mean"))


birest<- apply(bir.chain, c(1,2,3), mean)
dispest<- apply(disp.lik, c(1,2,3), mean)
par(mfrow=c(1,1))
plot(birest[,7,11] - italy.births.reg[,7,11], main = "Births: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", pch = 18)

dplot( ~ age | region, data = bir.chain, subarray = time == "2012",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation")

dplot( ~ age | region, data = disp.lik, subarray = time == "2012",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation")


dplot( ~ age | region, data = Counts(birest, dimscales = c(time = "Intervals")),
       subarray = time == "2012",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation",
       overlay = list(values = Counts(italy.births.reg[,,11]), col = "red", lwd= 2))


MCMC <-fetchMCMC(filename)

plot(MCMC$model.likelihood.count, sub= "Likelihood counts")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.hyper.age.scaleLevel, sub= "Hyper age scale level" )
plot(MCMC$model.prior.count.mean, sub= "Prior count mean" )
plot(MCMC$model.prior.count.sd, sub= "Prior count sd" )
plot(MCMC$model.hyper.age.scaleTrend, sub= "Hyper age scale trend")
plot(MCMC$model.hyper.age.scaleError, sub= "Hyper age scale error" )
plot(MCMC$model.hyper.year.damp, sub= "Hyper year damp" )
plot(MCMC$model.hyper.year.scaleError, sub= "Hyper year scale error" )

#------------------------------
# Births with exposure
#------------------------------
y <- italy.births.reg %>%
  Counts(dimscales = c(time = "Intervals"))

expose <- italy.popn.reg %>%
  Counts(dimscales = c(time = "Intervals")) %>%
  subarray(sex == "Female") %>%
  subarray(age >14 & age < 50)

filename <- tempfile()

estimateModel(Model(y ~ CMP(mean ~ age + region + time,
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 1),
                                                    scale = HalfT(scale = 0.1))),
                    age ~ DLM(damp = NULL),
                    jump = 0.02),
              y = y,
              exposure = expose,
              filename = filename,
              nBurnin = 300000,
              nSim = 300000,
              nChain = 3,
              nThin = 1000)

continueEstimation(filename,nBurnin = 30000,
                   nSim = 30000,
                   nThin = 750)

fetchSummary(filename)
showModel(filename)
listContents(filename)

bir.chain <- fetch(filename, where=c("model", "likelihood", "rate"))
disp.lik <- fetch(filename, where=c("model", "likelihood", "dispersion"))
disp.prior <- fetch(filename, where=c("model", "prior", "dispersion", "mean"))


birest<- apply(bir.chain, c(1,2,3), mean)
dispest <- apply(disp.lik, c(1,2,3), mean)

expose <- aperm(expose, c(2,1,3))
plot(birest[,,11]*expose[,,15] - italy.births.reg[,,11], main = "Births: estimation vs data 2015",
     xlab = "Regions", ylab = "Diff", pch = 18)

dplot( ~ age | region, data = bir.chain, subarray = time == "2012",
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation")

data <- Counts(birest[,,11]*expose[,,15])
dplot( ~ age | region, data = data,
       prob = c(0.025, 0.25, 0.5, 0.75, 0.975), scales = list(y = "free"), 
       main = "Births estimation",
       overlay = list(values = Counts(italy.births.reg[,,11]), col = "red", lwd= 2))


MCMC <-fetchMCMC(filename)

plot(MCMC$model.likelihood.rate, sub= "Likelihood rate")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.hyper.age.scaleLevel, sub= "Hyper age scale level" )
plot(MCMC$model.prior.rate.mean, sub= "Prior rate mean" )
plot(MCMC$model.prior.rate.sd, sub= "Prior rate sd" )
plot(MCMC$model.hyper.age.scaleTrend, sub= "Hyper age scale trend")
plot(MCMC$model.hyper.age.scaleError, sub= "Hyper age scale error" )

#------------------------------------
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

estimateModel(Model(y ~ CMP(mean ~ age + sex + year + age * sex + sex * year,
                            dispersion = Dispersion(mean = Norm(mean = 1, sd = 1),
                                                    scale = HalfT(scale = 1))),
                    age ~ DLM(damp = NULL),
                    jump = 0.02),
              y = y,
              exposure = expose,
              filename = filename,
              nBurnin = 5000,
              nSim = 5000,
              nChain = 4,
              nThin = 100)

continueEstimation(filename,nBurnin = 10000,
                   nSim = 10000,
                   nThin = 250)

fetchSummary(filename)

MCMC <-fetchMCMC(filename)
plot(MCMC$model.likelihood.rate, sub= "Likelihood rate")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.hyper.age.scaleLevel, sub= "Hyper age scale level" )
plot(MCMC$model.prior.rate.mean, sub= "Prior rate mean" )
plot(MCMC$model.prior.rate.sd, sub= "Prior rate sd" )
plot(MCMC$model.hyper.age.scaleTrend, sub= "Hyper age scale trend")
plot(MCMC$model.hyper.age.scaleError, sub= "Hyper age scale error" )
plot(MCMC$model.hyper.year.damp, sub= "Hyper year damp" )
plot(MCMC$model.hyper.year.scaleError, sub= "Hyper year scale error" )



#-------------------------------------------------------
# Italian data old update model
#-------------------------------------------------------

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

filename <- tempfile()
estimateModel(Model(y ~ CMP(mean ~ 1,
                            dispersion = Dispersion(mean = Norm(mean = 0, sd = 0.1),
                                                    scale = HalfT(scale = 0.1)),
                            useExpose = FALSE),
                    jump = 0.001),
              y = y,
              filename = filename,
              nBurnin = 20000,
              nSim = 20000,
              nChain = 3,
              nThin = 500)

fetchSummary(filename)


MCMC <-fetchMCMC(filename)
plot(MCMC$model.likelihood.count)
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main = "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main = "Prior dispersion sd" )
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
              nBurnin = 0,
              nSim = 200,
              nChain = 4,
              nThin = 1)
fetchSummary(filename)

MCMC <-fetchMCMC(filename)
plot(MCMC$model.likelihood.rate, sub= "Likelihood rate")
plot(MCMC$model.likelihood.dispersion, sub= "Likelihood dispersion" )
plot(MCMC$model.prior.dispersion.mean, main= "Prior dispersion mean")
plot(MCMC$model.prior.dispersion.sd, main= "Prior dispersion sd" )
plot(MCMC$model.hyper.age.scaleLevel, sub= "Hyper age scale level" )
plot(MCMC$model.prior.rate.mean, sub= "Prior rate mean" )
plot(MCMC$model.prior.rate.sd, sub= "Prior rate sd" )
plot(MCMC$model.hyper.age.scaleTrend, sub= "Hyper age scale trend")
plot(MCMC$model.hyper.age.scaleError, sub= "Hyper age scale error" )
plot(MCMC$model.hyper.year.damp, sub= "Hyper year damp" )
plot(MCMC$model.hyper.year.scaleError, sub= "Hyper year scale error" )


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
