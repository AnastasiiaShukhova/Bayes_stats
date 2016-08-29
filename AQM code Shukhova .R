#########################################################
#########          AQM FINAL PROJECT             ########
#########################################################
### Author: Anastasiia Shukhova
### Title: Is the road to hell paved with good intentions?
#            A study of democratic sanctions effectiveness


########### LOADING AND PREPARATION OF THE DATA #############


#set working directory
setwd("/Users/nastyashukhova/Dropbox/UniMannheim/PoliSci/SS_16/AQM/Final paper/Analysis and code")

##### WARNING: JAGS needs to be installed on the computer to package 'rjags' work. Mac users should install 3.4.0 version of JAGS
### you can do it here: 
# https://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/Mac%20OS%20X/JAGS-Mavericks-3.4.0.dmg/download


packages <- c("ggplot2", "foreign", "stringr", "pcse", "Matrix", "plm", "lme4", "rjags", "devtools", "gdata", "plyr","tidyr", 'ggthemes', "fda", "car",
              "ggmcmc", "data.table")
for (p in packages) {
  if (p %in% installed.packages()[,1]) require(p, character.only=T)
  else {
    install.packages(p)
    library(p, character.only=T)
  }
}

#load modules for rjags

library(rjags)
load.module("glm")
load.module("lecuyer")

# 
# ###  load the data
# 
# soestData  <- read.dta("/Users/nastyashukhova/Dropbox/UniMannheim/PoliSci/SS_16/AQM/Final\ paper/rereplication\ materials/Replication_Democratization copy.dta")
# SW  <- read.dta("/Users/nastyashukhova/Dropbox/UniMannheim/PoliSci/SS_16/AQM/Final\ paper/SoestWahman_2015_replication/SW_replicationJPR.dta")
# 
# 
# ###  add variables that might be usefull from the similar dataset 
# 
# soestData$westorgtie  <- SW$westorgtie 
# soestData$inflation  <- SW$wdiinflationgdp
# soestData$westTrade  <- SW$westtradelog
# soestData$ifhpol  <- SW$ifhpol

# save(soestData, file = "soestData.Rdata")
load("soestData.Rdata")
# delete countries that have never experienced a sanction otherwise varying coefficients for these countries are uninterpretable 

for (i in unique(soestData$cname)){
  if (all(is.na(soestData[soestData$cname == i,]$dm_sancgoal))){
    soestData  <- soestData[soestData$cname != i , ]
    next
  }
  if (sum(soestData[soestData$cname == i,]$dm_sancgoal, na.rm=T)==0){
    soestData  <- soestData[soestData$cname != i,]
    next
  }
}

# rename Pakistan
soestData[soestData$cname == "Pakistan (1972-)", ]$cname <- "Pakistan" 

#add an ID variable
soestData$numID  <- as.numeric(as.factor(soestData$cname))


# let's log GDP and trade because they are skewed and have very high variability and scale
# ggplot(data = soestData, aes(x =gdpconstant )) + geom_density()
soestData$logGdp  <- log(soestData$gdpconstant)

# do the same  with trade
soestData$logTrade  <- log(soestData$trade)    


#### Figure 1.  Number of introduced sanctions over the years (INTRODUCTION)  ###########

numberOfSanctions <- rep(NA, length(unique(soestData$year)) )
mean(numberOfSanctions)

# create a variable that corresponds to the number of sanctions per year
counter = 1
for (i in 1990: (1990+(length(unique(soestData$year))-1))){
  numberOfSanctions[counter] <- sum(soestData[soestData$year == i,]$dm_sancgoal, na.rm = T)
  counter = counter + 1
}

# get the total number of sanctions
sum(numberOfSanctions)

# plot the distribution 
dataPlot <- data.frame('year' = 1990:2010, numberOfSanctions)
plot_numSanctions <- ggplot(data =dataPlot, aes(x = 1990:2010, y = numberOfSanctions)) + geom_point() + geom_line() +
  theme_tufte() + xlab("Year") + ylab("Number of Sanctions per Year") + scale_x_continuous(breaks = seq(1990, 2010, by = 2)) +
  scale_y_continuous(breaks = seq(8, 17, by = 1))

#save the plot
ggsave("plot_numSanctions.pdf", width = 7, height = 5)




############## DATA ANALYSIS ################
# exploit hierarchical centering 
soestData.dt <- as.data.table(soestData)

# function that centers variables by countries
centerWithinCountries <- function(data,  variable, group){
  data.dt = as.data.table(data)
  meanByGroup <- data.dt[,.(meanByGroup=mean(variable, na.rm = T)), by = group]
  centeredVar <- c()
  for (numb in unique(group)){
    for (i in variable[group == numb]){
      centeredVar <- rbind(centeredVar, i - as.numeric(meanByGroup[,.(meanByGroup)][numb]) )
    }
  }
  return(centeredVar)
}


# hierarchical centering 

soestData$centeredPop <- centerWithinCountries(data =  soestData, variable = soestData$logpwt_pop2, group= soestData$numID)
soestData$centeredGDP <- centerWithinCountries(data =  soestData, variable = soestData$logGdp, group= soestData$numID)
soestData$centeredTrade <- centerWithinCountries(data =  soestData, variable = soestData$logTrade, group= soestData$numID)
soestData$centeredWestTrade <- centerWithinCountries(data =  soestData, variable = soestData$westTrade, group= soestData$numID)
soestData$centeredOil <- centerWithinCountries(data =  soestData, variable = soestData$oilmil, group= soestData$numID)

#### Model where I account for the time dependencies using the Jackman's (2009) method #####
cat('model{
    # MODEL
    for(i in 1:N){
    y[i] ~ dnorm(y.hat[i], tau)
    y.hat[i] <- inprod(B[country[i], ], X[i,]) + inprod(b2, X.0[i,]) 
    for (j in 1:K){
    X[i,j]  ~ dnorm(0,  .0001)
    }
    for (k in 1:K.0){
    X.0[i,k]  ~ dnorm(0,  .0001)
    }
    }
    tau ~ dgamma(.1,.1)
    
    
    # PRIORS
    
    # set priors for varying intercepts and varying slopes allowing for the correlation between them 
    for (j in 1:J){ 
    for (k in 1:K){
    B[j,k] <- xi[k]*B.raw[j,k] 
    }
    B.raw[j,1:K] ~ dmnorm(mu.raw[], Tau.B.raw[,]) 
    }
    for (k in 1:K){
    mu[k] <- xi[k]*mu.raw[k] 
    mu.raw[k] ~ dnorm (0, .0001) 
    xi[k] ~ dunif (0, 100)
    }
    
    # correlation matrix of slopes with prior from Wishart distribution
    Tau.B.raw[1:K,1:K] ~ dwish (W[,], df)
    df <- K+1
    Sigma.B.raw[1:K,1:K] <- inverse(Tau.B.raw[,]) 
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.b[k,k.prime] <- Sigma.B.raw[k,k.prime]/sqrt(Sigma.B.raw[k,k]*Sigma.B.raw[k.prime,k.prime]) 
    }
    sigma.B[k] <- abs(xi[k])*sqrt(Sigma.B.raw[k,k]) 
    }
    
    for (k in 1:K.0){
    b2[k] ~ dnorm (0, taub2)
    }
    taub2 ~ dgamma(.1,.01)
    }',file={t.jackman <- tempfile()})


# set data for the model
X  <- cbind("beta0" = 1, "yearDif" = soestData$year-mean(soestData$year), "laggedSanction" = stats::lag(soestData$dm_sancgoal))

X0  <- data.frame( "Population" = stats::lag(soestData$centeredPop), "log GDP" = stats::lag(soestData$centeredGDP), "logTrade" = stats::lag(soestData$centeredTrade), 
                   "Civil War" = stats::lag(soestData$civilwar), "West Trade" = stats::lag(soestData$centeredWestTrade),
                   "Protest" = stats::lag(soestData$protest), "Oil" = stats::lag(soestData$centeredOil))

# parameter for the Wishart distribution
W <- diag (3)

# create data for the list for JAGS 
data.jackman <- list("y"=soestData$ifhpol,
                     "X"=X,
                     "country"=as.numeric(as.factor(soestData$cname)),
                     "N"=nrow(soestData),
                     "J"=length(unique(soestData$cname)),
                     "X.0" = X0,
                     "K.0" = ncol(X0),
                     "K" = ncol(X),
                     "W" = W
)


# set inits
random.no <- round(runif(3,0,1000),0)
inits <- list( 
  list(".RNG.name"="lecuyer::RngStream", ".RNG.seed"=random.no[1]),
  list(".RNG.name"="lecuyer::RngStream", ".RNG.seed"=random.no[2]),
  list(".RNG.name"="lecuyer::RngStream", ".RNG.seed"=random.no[3])
)


# set parameters that we want to save 
parameters.to.save <-c ("B", "b2", "mu", "sigma.B","tau", "rho.b","taub2")


# MODEL
jags.jackman  <- jags.model(t.jackman, # temporary file
                            data=data.jackman ,
                            inits=inits,
                            n.chain=3
)


# SAMPLING

####   W-A-R-N-I-N-G!!!! on my computer it takes 30-40 minutes to calculate this 
# please load already calculated MCMC object to save your time

mcmc.results.jackman <- coda.samples(model=jags.jackman,
                                     variable.names=parameters.to.save,
                                     n.iter=40000,
                                     n.adapt = 20000,
                                     thin = 40)


# save(mcmc.results.jackman, file = "mcmc.results.jackman.Rdata")

#load calculated model
load("mcmc.results.jackman.Rdata")



######   Analysis of the results   #######

#add id variable
soestData$numID  <- as.numeric(as.factor(soestData$cname))


# get the coeffients and quantiles
res.jackman <-  summary(mcmc.results.jackman)$quantiles
res.jackman.fixed <- summary(mcmc.results.jackman)$statistics


# get number of countries with significan effects
numOfCountries  <- length(unique(soestData$cname))
significant.jackman  <- c()
for (i in 1:numOfCountries){
  if ((res.jackman[((numOfCountries*2) + i), 1] < 0) & (res.jackman[((numOfCountries*2) + i), 5] < 0)){
    significant.jackman  <- rbind(significant.jackman, i)
  }
}
significantCountries.jackman <- rep(NA, length(significant.jackman))



counter = 1 
for (i in significant.jackman){
  significantCountries.jackman[counter]  <- unique(soestData[soestData$numID == i, ]$cname)
  counter <- counter + 1 
}



#### Presentation of the results  ######

# save the mcmc object as an ggmcmc object (for the ggmcmc package)
# S <- ggs(mcmc.results.jackman, burnin = 20000)
# # get the .pdf file with different convergence plots WARNING: takes some time to compile 
# ggmcmc(S, family = "B")
# ggmcmc(S, family = "b", file = "controls.pdf")


# get nice plots for the coefficients 
library(data.table)

soestData.dt <- as.data.table(soestData)
demMean <- soestData.dt[,.(meanDem=median(ifhpol, na.rm = T)), by = cname]

L.sanction.coef <- data.frame(
  Parameter=paste("B[", unique(soestData$numID), ",3]", sep=""),
  Label=unique(soestData$cname),
  MeanDem = demMean[,.(meanDem)]
)

suppressWarnings(S.full <- ggs(mcmc.results.jackman[,c(79:117)], 
              par_labels=L.sanction.coef, family="^B"))

results <- ggs_caterpillar(S.full) + theme_bw() + 
  aes(color=L.sanction.coef$meanDem) + 
  scale_color_continuous(guide_legend(title= "Median Level \nof Democracy")) +
  ylab("Effect of Democratic Sanctions") +theme_tufte(ticks  = F) +
  geom_vline(data=NULL, aes(xintercept=0), color="darkblue", size=0.25)
  
quartz()
results

# save plot 
ggsave("results.pdf", height = 8, width = 10)

### Figure 3. The posterior distribution of the effect of the democratic sanctions on the level of democracy in Belarus

EffectBelarus <- data.frame(
  Parameter="B[2,3]",
  Label= " ")


suppressMessages(suppressWarnings(
  plotBelarus <-  ggs_density(ggs(mcmc.results.jackman, par_labels=EffectBelarus, family="B\\[2,3\\]", burnin = 20000 )) + theme_tufte() + 
  scale_color_manual(values = c("#CCCCCC", "#999999", "#666666"), guide = guide_legend(title= "Chain")) +
  scale_fill_manual(values = c("#CCCCCC", "#999999", "#666666"), guide = guide_legend(title= "Chain")) + xlab("Change in the level of democracy after the \n introduction of democratic sanctions in Belarus"))) 

plotBelarus <- plotBelarus + geom_vline(xintercept = c(as.numeric(res.jackman[79+1, c(1, 3, 5)])), size = 0.25, linetype = "dashed", color = "#333333") 

quartz()
plotBelarus
# save plot 
ggsave("plotBelarus.pdf", plot = plotBelarus, height = 2, width = 4)




### Figure 4. The posterior distribution of the effect of the democratic sanctions on the level of democracy in Haiti

EffectHaiti <- data.frame(
  Parameter="B[17,3]",
  Label= " ")

suppressMessages(suppressWarnings(
plotHaiti <-  ggs_density(ggs(mcmc.results.jackman, par_labels=EffectHaiti, family="\\[17,3\\]", burnin = 20000 )) + theme_tufte() + 
  scale_color_manual(values = c("#CCCCCC", "#999999", "#666666"), guide = guide_legend(title= "Chain")) +
  scale_fill_manual(values = c("#CCCCCC", "#999999", "#666666"), guide = guide_legend(title= "Chain")) + xlab("Change in the level of democracy after the \n introduction of democratic sanctions in Haiti") ))

plotHaiti <- plotHaiti + geom_vline(xintercept = c(as.numeric(res.jackman[79+16, c(1, 3, 5)])), size = 0.25, linetype = "dashed", color = "#333333") 

quartz()
plotHaiti

# save plot 
ggsave("plotHaiti.pdf", plot = plotHaiti, height = 2, width = 4)




###################    APPENDICES         #################

####### Appendix A. #####
# LaTeX output for the descriptive statistics
library(stargazer)
stargazer(soestData[,c("ifhpol", "dm_sancgoal", "westTrade", "logGdp", "logTrade", "oilmil")])

####### Appendix B. #####

#### plot level of democracy over the years by country

plot_countries_over_time <- ggplot(data =soestData, aes(x = year, y = ifhpol)) +
  geom_point() + geom_line() + facet_wrap(~cname, 5) + theme_bw() + ylab("Level of Democracy") + xlab("Year")

quartz()
plot_countries_over_time

# save plot 
ggsave(filename =" dem_over_time.pdf", plot = plot_countries_over_time, width = 25, height = 15)


####### Appendix D. #####

rhats <- gelman.diag(mcmc.results.jackman[,c(79:117)])

suppressMessages(rhats.hist <- ggplot(data = NULL, aes(x = rhats$psrf[,1])) + 
                   geom_histogram() + xlab("Gelman's R.hat") + theme_tufte() )
suppressMessages(ggsave("rhats.hist.pdf", plot = rhats.hist, height = 4, width = 4))


####### Appendix E. #####

# we can use another method for accounting for the time trend
# # add Voice and accointability 
library("readstata13", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
vdem <- read.dta13("VDem.dta")

# rename variables so that they have the same namens of countries and could be easily merged 

vdem[vdem$country_name=="Congo_Democratic Republic of",]$country_name <- "Congo, Democratic Republic"
vdem[vdem$country_name=="Korea_North",]$country_name <- "Korea, North"

soestData.vdem  <- merge(soestData, vdem, by.x = c("cname", "year"), by.y = c("country_name", "year"))
# length(unique(soestData$cname))



#### Accounting for the time dependencies using the Jackman's method #####
cat('model{
    # MODEL
    for(i in 1:N){
    y[i] ~ dnorm(y.hat[i], tau)
    y.hat[i] <- inprod(B[country[i], ], X[i,]) + inprod(b2, X.0[i,]) 
    for (j in 1:K){
    X[i,j]  ~ dnorm(0,  .0001)
    }
    for (k in 1:K.0){
    X.0[i,k]  ~ dnorm(0,  .0001)
    }
    }
    tau ~ dgamma(.1,.1)
    
    # PRIORS
    
    # set priors for varying intercepts and varying slopes allowing for the correlation between them 
    for (j in 1:J){ 
    for (k in 1:K){
    B[j,k] <- xi[k]*B.raw[j,k] 
    }
    B.raw[j,1:K] ~ dmnorm(mu.raw[], Tau.B.raw[,]) 
    }
    for (k in 1:K){
    mu[k] <- xi[k]*mu.raw[k] 
    mu.raw[k] ~ dnorm (0, .0001) 
    xi[k] ~ dunif (0, 100)
    }
    
    # correlation matrix of slopes with prior from Wishart distribution
    Tau.B.raw[1:K,1:K] ~ dwish (W[,], df)
    df <- K+1
    Sigma.B.raw[1:K,1:K] <- inverse(Tau.B.raw[,]) 
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.b[k,k.prime] <- Sigma.B.raw[k,k.prime]/sqrt(Sigma.B.raw[k,k]*Sigma.B.raw[k.prime,k.prime]) 
    }
    sigma.B[k] <- abs(xi[k])*sqrt(Sigma.B.raw[k,k]) 
    }
    
    for (k in 1:K.0){
    b2[k] ~ dnorm (0, taub2)
    }
    taub2 ~ dgamma(.1,.01)
    }',file={t.jackman <- tempfile()})

X  <- cbind("beta0" = 1, "yearDif" = soestData.vdem$year-mean(soestData.vdem$year), "laggedSanction" = stats::lag(soestData.vdem$dm_sancgoal))

X0  <- data.frame( "Population" = stats::lag(soestData.vdem$centeredPop), "log GDP" = stats::lag(soestData.vdem$centeredGDP), "logTrade" = stats::lag(soestData.vdem$centeredTrade), 
                   "Civil War" = stats::lag(soestData.vdem$civilwar), "West Trade" = stats::lag(soestData.vdem$centeredWestTrade),
                   "Protest" = stats::lag(soestData.vdem$protest), "Oil" = stats::lag(soestData.vdem$oilProd))


W <- diag (3)

data.jackman.vdem <- list("y"=soestData.vdem$v2x_polyarchy,
                          "X"=X,
                          "country"=as.numeric(as.factor(soestData.vdem$cname)),
                          "N"=nrow(soestData.vdem),
                          "J"=length(unique(soestData.vdem$cname)),
                          "X.0" = X0,
                          "K.0" = ncol(X0),
                          "K" = ncol(X),
                          "W" = W
)


# set inits
random.no <- round(runif(3,0,1000),0)
inits <- list( 
  list(".RNG.name"="lecuyer::RngStream", ".RNG.seed"=random.no[1]),
  list(".RNG.name"="lecuyer::RngStream", ".RNG.seed"=random.no[2]),
  list(".RNG.name"="lecuyer::RngStream", ".RNG.seed"=random.no[3])
)

# set parameters that we want to save 
parameters.to.save <-c ("B", "b2", "mu", "sigma.B","tau", "rho.b","taub2")

# MODEL
jags.jackman.vdem  <- jags.model(t.jackman, # temporary file
                                 data=data.jackman.vdem ,
                                 inits=inits,
                                 n.chain=3
)


####   W-A-R-N-I-N-G!!!! on my computer it takes 30-40 minutes to calculate this 
# please load already calculated MCMC object to save your time

# SAMPLING
mcmc.results.vdem <- coda.samples(model=jags.jackman.vdem,
                                  variable.names=parameters.to.save,
                                  n.iter=15000,
                                  n.adapt = 10000,
                                  thin = 30)


# save(mcmc.results.vdem, file = "mcmc.results.vdem.Rdata")
load("mcmc.results.vdem.Rdata")
result.vdem  <-  summary(mcmc.results.vdem)$quantiles


# get nice plots for the coefficients 
soestData.vdem$numID  <- as.numeric(as.factor(soestData.vdem$cname))
data.frame(unique(soestData.vdem$numID), unique(soestData.vdem$cname))

L.sanction.coef <- data.frame(
  Parameter=paste("B[", unique(soestData.vdem$numID), ",3]", sep=""),
  Label=unique(soestData.vdem$cname)
)
mcmc.results.vdem.clean <- mcmc.results.vdem[,((36*2)+1):(36*3)]
S.full <- ggs(mcmc.results.vdem.clean[,-c(18, 9 , 6)], 
              par_labels=L.sanction.coef, family="B")


results.check <- ggs_caterpillar(S.full) + theme_bw() + 
  ylab("Effect of Democratic Sanctions") +theme_tufte(ticks  = F) + geom_vline(data=NULL, aes(xintercept=0), color="darkblue", size=0.25)

quartz()
results.check
ggsave("results.check.pdf", height = 8, width = 10)

