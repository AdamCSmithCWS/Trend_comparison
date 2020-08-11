### example of how to compare estimates of population trends while accounting for the precision of those trend estimates

library(googledrive)
library(tidyverse)
library(jagsUI)
library(ggmcmc)


# Setup downloading the 2018 BBS trends from my Google Drive --------------
### this commented out section doesn't need to be re-run
# tdwn <- drive_download(file = "All 2018 BBS trends.csv",
#                        path = "data/All_2018_BBS_trends.csv")
# 
# 
# trends <- read.csv("data/All_2018_BBS_trends.csv",stringsAsFactors = F)
#

# ### selecting out the continental long-term trends and the trend-relevant columns
tr1 = trends %>%
  filter(Region_type == "national", Trend_Time == "Long-term") %>%
  select(species,Trend, Trend_Q0.025, Trend_Q0.975,Number_of_Routes,Region)

write.csv(tr1,"data/national_longterm_trends.csv")





# reading in continental BBS trends ---------------------------------------

tr = read.csv("data/national_longterm_trends.csv")

### funciton to transform %-change trends to log-scale geometric mean rates of change
log_trans <- function(x){
  log((x/100)+1)
}

## transforming the 95% CIs into an estimate of precision (1/variance)
##pivoting wider to separate US and Canadian estimates
##filtering to species that have estimates for both
   trlong <- tr %>% mutate(betahat = log_trans(Trend),
           log_lci = log_trans(Trend_Q0.025),
           log_uci = log_trans(Trend_Q0.975)) %>%
    mutate(varhat = ((log_uci-log_lci)/(1.96*2))^2) %>% 
     pivot_wider(id_cols = species,names_from = Region, values_from = c(betahat,varhat,Number_of_Routes) ) %>% 
     filter(.,complete.cases(betahat_US,betahat_CA))

 
   nspecies = nrow(trlong)
   nroutes = as.matrix(trlong[,c("Number_of_Routes_US","Number_of_Routes_CA")])
   betahat = as.matrix(trlong[,c("betahat_US","betahat_CA")])
   varhat = as.matrix(trlong[,c("varhat_US","varhat_CA")])
   
   jags_data <- list(nspecies = nspecies,
                     nroutes = nroutes,
                     betahat = betahat,
                     varhat = varhat)
# model -------------------------------------------------------------------

 modl <- "

### model for comparing trend estimates from two (indexed by e) different analyses (or regions), for species (indexed by s)
### for this comparison, e==1 are estimates from the US and e==2 from Canada
### Model is modified from the grassland bird model from Link and Barker 2010, pg 90, also described in Sauer and Link 2002, Ecology 83:1743-1751
### in essence, this is a paired t-test style comparison, that accounts for the imprecision in each s*e trend estimate.

### input data compiled in R, consist of nspecies (the number of species), as well as 3 matrices: varhat, betahat, and n, each of which has nspecies rows and 2 columns
## varhat = estimates of the log-scale variances of trends
## betahat = estimates of the log-scale trends
## nroutes = sample size of the trends - number of routes used to generate the trend)
## nspecies = number of species

model{

for (e in 1:2) {

for(s in 1:nspecies) {

	varhat[s,e] ~ dgamma(p[s,e],lam[s,e]) #varhat = data = variance of trend estimates
	p[s,e] <- nroutes[s,e] / 2
	lam[s,e] <- p[s,e] * tau.betahat[s,e]
	tau.betahat[s,e] ~ dgamma(0.001,0.001)
	sd.betahat[s,e] <- 1/sqrt(tau.betahat[s,e])
	### above is prior for the estimated species level precision (tau-betahat), accounting fo rhte number of routes included
	
	betahat[s,e] ~ dnorm(beta[s,e],tau.betahat[s,e]) #betahat = data = trend estimates
	beta[s,e] ~ dnorm(0.0,1)
	#### beta is the prior for the estimated species level trend, accounting for the precision
	
	 pos[s,e] <- step(beta[s,e]) #tracks whether the trend is positive
	 neg[s,e] <- 1-pos[s,e] #tracks opposite
	

} # end of first s-species loop
  numpos[e] <- sum(pos[,e]) # if pos is tracked, then this estimates the number of species with positive trends
  numneg[e] <- sum(neg[,e]) # if neg is tracked, same

	} #end of e loop (indexing two models being compared)

### species level differences between US and Canadian trends (US - Canadian)
### if dif[s] is positive = US trend is more positive than Canadian trend
for(s in 1:nspecies) {
dif[s] <- beta[s,1]-beta[s,2] # dif is a vector of the species-specific trend differences after accounting for the imprecision of each estimate's trend and the group (survey/monitoring program) structure
difvar[s] <- sd.betahat[s,1]-sd.betahat[s,2] # difvar is a vector of the species-specific differences in sd of the trends
} # end of second s-species loop


m.dif <- mean(dif[])
m.difvar <- mean(difvar[])

dif.numpos <- numpos[1]-numpos[2]
dif.numneg <- numneg[1]-numneg[2]


} # end of model

   "  
   
   trend_comp = "trend_comparison.R"
cat(modl,file = trend_comp)   


params <- c("m.dif",
            "m.difvar",
            #"dif",
            #"beta",
            #"difvar",
            "dif.numneg")


burnInSteps = 2000            # Number of steps to "burn-in" the samplers. this is sufficient for testing, but you'll want to increase this
nChains = 3                   # Number of chains to run.
numSavedSteps=1000         # Total number of steps in each chain to save. this is sufficient for testing, but you'll want to increase this
thinSteps=10                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.



out = jagsUI(data = jags_data,
             parameters.to.save = params,
             n.chains = 3,
             n.burnin = burnInSteps,
             n.thin = thinSteps,
             n.iter = nIter,
             parallel = T,
             model.file = trend_comp)



summr = out$summary #table showing a summary of the posterior distribution and some basic convergence stats for all the monitored parameters

out_ggs = ggs(out$samples)
ggmcmc(out_ggs,file = "convergence summaries.pdf", paparam_page = 8)



